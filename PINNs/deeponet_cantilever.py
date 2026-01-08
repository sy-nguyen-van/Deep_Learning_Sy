"""
DeepONet for 2D Cantilever Beam Linear Elasticity

This implementation uses Deep Operator Networks to learn the solution operator
mapping from load distributions to displacement fields for a cantilever beam.

Problem Setup:
- Cantilever beam: Length=4m, Height=1m
- Material: E=210 GPa, ν=0.3 (steel-like)
- Boundary conditions: Fixed on left edge, free on right edge
- Variable loads on top surface
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple, List, Optional
import os


# ============================================================================
# DeepONet Architecture
# ============================================================================

class BranchNetwork(nn.Module):
    """Branch network processes the input function (load distribution)"""
    
    def __init__(self, sensor_count: int, hidden_dims: List[int], output_dim: int):
        super(BranchNetwork, self).__init__()
        
        layers = []
        input_dim = sensor_count
        
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, hidden_dim))
            layers.append(nn.Tanh())
            input_dim = hidden_dim
        
        layers.append(nn.Linear(input_dim, output_dim))
        
        self.network = nn.Sequential(*layers)
    
    def forward(self, u):
        """
        Args:
            u: Input function values at sensor locations [batch_size, sensor_count]
        Returns:
            Branch output [batch_size, output_dim]
        """
        return self.network(u)


class TrunkNetwork(nn.Module):
    """Trunk network processes the spatial coordinates"""
    
    def __init__(self, input_dim: int, hidden_dims: List[int], output_dim: int):
        super(TrunkNetwork, self).__init__()
        
        layers = []
        
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, hidden_dim))
            layers.append(nn.Tanh())
            input_dim = hidden_dim
        
        layers.append(nn.Linear(input_dim, output_dim))
        
        self.network = nn.Sequential(*layers)
    
    def forward(self, y):
        """
        Args:
            y: Spatial coordinates [batch_size, n_points, input_dim]
        Returns:
            Trunk output [batch_size, n_points, output_dim]
        """
        return self.network(y)


class DeepONet(nn.Module):
    """
    Deep Operator Network for learning solution operators
    
    The network learns G: u → s where:
    - u is the input function (load distribution)
    - s is the output function (displacement field)
    """
    
    def __init__(
        self,
        sensor_count: int,
        branch_hidden_dims: List[int],
        trunk_hidden_dims: List[int],
        latent_dim: int,
        n_outputs: int = 2  # u_x and u_y
    ):
        super(DeepONet, self).__init__()
        
        self.n_outputs = n_outputs
        self.latent_dim = latent_dim
        
        # Create separate DeepONets for each output component
        self.branch_nets = nn.ModuleList([
            BranchNetwork(sensor_count, branch_hidden_dims, latent_dim)
            for _ in range(n_outputs)
        ])
        
        self.trunk_nets = nn.ModuleList([
            TrunkNetwork(2, trunk_hidden_dims, latent_dim)  # 2D coordinates
            for _ in range(n_outputs)
        ])
        
        # Bias terms
        self.biases = nn.ParameterList([
            nn.Parameter(torch.zeros(1))
            for _ in range(n_outputs)
        ])
    
    def forward(self, u, y):
        """
        Args:
            u: Input function at sensors [batch_size, sensor_count]
            y: Query points [batch_size, n_points, 2]
        Returns:
            Output function values [batch_size, n_points, n_outputs]
        """
        batch_size, n_points, _ = y.shape
        outputs = []
        
        for i in range(self.n_outputs):
            # Branch network output
            branch_out = self.branch_nets[i](u)  # [batch_size, latent_dim]
            
            # Trunk network output
            y_flat = y.reshape(-1, 2)  # [batch_size * n_points, 2]
            trunk_out = self.trunk_nets[i](y_flat)  # [batch_size * n_points, latent_dim]
            trunk_out = trunk_out.reshape(batch_size, n_points, self.latent_dim)
            
            # Inner product
            # branch_out: [batch_size, latent_dim]
            # trunk_out: [batch_size, n_points, latent_dim]
            output = torch.sum(
                branch_out.unsqueeze(1) * trunk_out,
                dim=2
            ) + self.biases[i]  # [batch_size, n_points]
            
            outputs.append(output)
        
        # Stack outputs
        result = torch.stack(outputs, dim=2)  # [batch_size, n_points, n_outputs]
        return result


# ============================================================================
# Data Generation
# ============================================================================

class CantileverBeamDataGenerator:
    """
    Generate training data for cantilever beam using analytical solutions
    or simplified FEM approach
    """
    
    def __init__(
        self,
        length: float = 4.0,
        height: float = 1.0,
        E: float = 210e9,  # Young's modulus (Pa)
        nu: float = 0.3,   # Poisson's ratio
        n_sensors: int = 50
    ):
        self.L = length
        self.H = height
        self.E = E
        self.nu = nu
        self.n_sensors = n_sensors
        
        # Sensor locations (along the top surface)
        self.sensor_x = np.linspace(0, self.L, n_sensors)
        self.sensor_y = np.ones(n_sensors) * self.H
        
        # Material properties
        self.D = self.E / (1 - self.nu**2)  # Plane stress
    
    def generate_random_load(self) -> np.ndarray:
        """Generate random load distribution on top surface"""
        # Random combination of basis functions
        n_modes = np.random.randint(1, 5)
        load = np.zeros(self.n_sensors)
        
        for _ in range(n_modes):
            mode_type = np.random.choice(['constant', 'linear', 'sine', 'point'])
            amplitude = np.random.uniform(-1e6, -1e4)  # Downward loads
            
            if mode_type == 'constant':
                load += amplitude * np.ones(self.n_sensors)
            elif mode_type == 'linear':
                slope = np.random.uniform(-1, 1)
                load += amplitude * (1 + slope * self.sensor_x / self.L)
            elif mode_type == 'sine':
                freq = np.random.randint(1, 4)
                load += amplitude * np.sin(freq * np.pi * self.sensor_x / self.L)
            elif mode_type == 'point':
                center = np.random.uniform(0.2 * self.L, 0.8 * self.L)
                width = np.random.uniform(0.1 * self.L, 0.3 * self.L)
                load += amplitude * np.exp(-((self.sensor_x - center) / width)**2)
        
        return load
    
    def solve_beam_deflection(self, load: np.ndarray, x_eval: np.ndarray, y_eval: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simplified analytical solution for beam deflection
        Using Euler-Bernoulli beam theory for vertical displacement
        and simplified horizontal displacement
        
        Args:
            load: Load distribution at sensor points [n_sensors]
            x_eval: x-coordinates for evaluation [n_points]
            y_eval: y-coordinates for evaluation [n_points]
        
        Returns:
            u_x, u_y: Displacement components [n_points]
        """
        # Moment of inertia
        I = (self.H**3) / 12
        
        # Integrate load to get shear and moment
        # Approximate load as piecewise constant
        dx = self.L / (self.n_sensors - 1)
        
        u_x = np.zeros_like(x_eval)
        u_y = np.zeros_like(x_eval)
        
        for i, (x, y) in enumerate(zip(x_eval, y_eval)):
            # Vertical displacement using beam theory
            # v(x) = integral of M(x) / EI
            
            # Calculate moment at position x
            M = 0
            V = 0
            for j, (x_load, q) in enumerate(zip(self.sensor_x, load)):
                if x_load > x:
                    # Load to the right contributes to moment
                    V += q * dx
                    M += q * dx * (x_load - x)
            
            # Deflection (simplified)
            # Using v''(x) = M(x) / EI and integrating
            u_y[i] = M * x**2 / (2 * self.E * I)
            
            # Add contribution from shear
            u_y[i] += V * x**3 / (6 * self.E * I)
            
            # Horizontal displacement (simplified, much smaller)
            # Due to Poisson effect and bending
            u_x[i] = -self.nu * (y - self.H/2) * M * x / (self.E * I)
        
        return u_x, u_y
    
    def generate_sample(self, n_points: int = 100) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate one training sample
        
        Returns:
            load: Load at sensor points [n_sensors]
            coords: Evaluation coordinates [n_points, 2]
            displacements: Displacement field [n_points, 2]
        """
        # Generate random load
        load = self.generate_random_load()
        
        # Random evaluation points in the domain
        x_eval = np.random.uniform(0, self.L, n_points)
        y_eval = np.random.uniform(0, self.H, n_points)
        
        # Solve for displacements
        u_x, u_y = self.solve_beam_deflection(load, x_eval, y_eval)
        
        coords = np.stack([x_eval, y_eval], axis=1)
        displacements = np.stack([u_x, u_y], axis=1)
        
        return load, coords, displacements
    
    def generate_dataset(self, n_samples: int, n_points_per_sample: int = 100):
        """Generate dataset for training"""
        loads = []
        coords = []
        displacements = []
        
        for _ in range(n_samples):
            load, coord, disp = self.generate_sample(n_points_per_sample)
            loads.append(load)
            coords.append(coord)
            displacements.append(disp)
        
        return (
            np.array(loads),
            np.array(coords),
            np.array(displacements)
        )


class CantileverDataset(Dataset):
    """PyTorch Dataset for cantilever beam data"""
    
    def __init__(self, loads, coords, displacements):
        self.loads = torch.FloatTensor(loads)
        self.coords = torch.FloatTensor(coords)
        self.displacements = torch.FloatTensor(displacements)
    
    def __len__(self):
        return len(self.loads)
    
    def __getitem__(self, idx):
        return self.loads[idx], self.coords[idx], self.displacements[idx]


# ============================================================================
# Training
# ============================================================================

class DeepONetTrainer:
    """Trainer for DeepONet"""
    
    def __init__(
        self,
        model: DeepONet,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu'
    ):
        self.model = model.to(device)
        self.device = device
        self.history = {'train_loss': [], 'val_loss': []}
    
    def train_epoch(self, dataloader, optimizer, criterion):
        """Train for one epoch"""
        self.model.train()
        total_loss = 0
        
        for loads, coords, displacements in dataloader:
            loads = loads.to(self.device)
            coords = coords.to(self.device)
            displacements = displacements.to(self.device)
            
            optimizer.zero_grad()
            
            # Forward pass
            predictions = self.model(loads, coords)
            
            # Compute loss
            loss = criterion(predictions, displacements)
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
        
        return total_loss / len(dataloader)
    
    def validate(self, dataloader, criterion):
        """Validate the model"""
        self.model.eval()
        total_loss = 0
        
        with torch.no_grad():
            for loads, coords, displacements in dataloader:
                loads = loads.to(self.device)
                coords = coords.to(self.device)
                displacements = displacements.to(self.device)
                
                predictions = self.model(loads, coords)
                loss = criterion(predictions, displacements)
                
                total_loss += loss.item()
        
        return total_loss / len(dataloader)
    
    def train(
        self,
        train_loader,
        val_loader,
        n_epochs: int,
        learning_rate: float = 1e-3,
        scheduler_step: int = 100,
        scheduler_gamma: float = 0.9
    ):
        """Full training loop"""
        optimizer = optim.Adam(self.model.parameters(), lr=learning_rate)
        scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=scheduler_step, gamma=scheduler_gamma)
        criterion = nn.MSELoss()
        
        print(f"Training on device: {self.device}")
        print(f"Number of parameters: {sum(p.numel() for p in self.model.parameters())}")
        
        for epoch in range(n_epochs):
            train_loss = self.train_epoch(train_loader, optimizer, criterion)
            val_loss = self.validate(val_loader, criterion)
            
            self.history['train_loss'].append(train_loss)
            self.history['val_loss'].append(val_loss)
            
            scheduler.step()
            
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch+1}/{n_epochs} - Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}")
        
        print("Training completed!")
    
    def save_model(self, path: str):
        """Save model checkpoint"""
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'history': self.history
        }, path)
        print(f"Model saved to {path}")
    
    def load_model(self, path: str):
        """Load model checkpoint"""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.history = checkpoint['history']
        print(f"Model loaded from {path}")


# ============================================================================
# Visualization
# ============================================================================

class Visualizer:
    """Visualization utilities for cantilever beam results"""
    
    def __init__(self, L: float = 4.0, H: float = 1.0):
        self.L = L
        self.H = H
    
    def plot_displacement_field(
        self,
        model: DeepONet,
        load: np.ndarray,
        device: str = 'cpu',
        resolution: int = 50,
        save_path: Optional[str] = None
    ):
        """Plot displacement field for a given load"""
        model.eval()
        
        # Create grid
        x = np.linspace(0, self.L, resolution)
        y = np.linspace(0, self.H, resolution)
        X, Y = np.meshgrid(x, y)
        
        coords = np.stack([X.flatten(), Y.flatten()], axis=1)
        coords_tensor = torch.FloatTensor(coords).unsqueeze(0).to(device)
        load_tensor = torch.FloatTensor(load).unsqueeze(0).to(device)
        
        # Predict
        with torch.no_grad():
            disp = model(load_tensor, coords_tensor)
            disp = disp.cpu().numpy()[0]
        
        u_x = disp[:, 0].reshape(resolution, resolution)
        u_y = disp[:, 1].reshape(resolution, resolution)
        u_mag = np.sqrt(u_x**2 + u_y**2)
        
        # Plot
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        # u_x
        im0 = axes[0].contourf(X, Y, u_x, levels=20, cmap='RdBu_r')
        axes[0].set_title('Horizontal Displacement $u_x$')
        axes[0].set_xlabel('x (m)')
        axes[0].set_ylabel('y (m)')
        axes[0].set_aspect('equal')
        plt.colorbar(im0, ax=axes[0], label='Displacement (m)')
        
        # u_y
        im1 = axes[1].contourf(X, Y, u_y, levels=20, cmap='RdBu_r')
        axes[1].set_title('Vertical Displacement $u_y$')
        axes[1].set_xlabel('x (m)')
        axes[1].set_ylabel('y (m)')
        axes[1].set_aspect('equal')
        plt.colorbar(im1, ax=axes[1], label='Displacement (m)')
        
        # Magnitude
        im2 = axes[2].contourf(X, Y, u_mag, levels=20, cmap='viridis')
        axes[2].set_title('Displacement Magnitude')
        axes[2].set_xlabel('x (m)')
        axes[2].set_ylabel('y (m)')
        axes[2].set_aspect('equal')
        plt.colorbar(im2, ax=axes[2], label='Displacement (m)')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        
        plt.show()
    
    def plot_load_distribution(self, load: np.ndarray, sensor_x: np.ndarray, save_path: Optional[str] = None):
        """Plot the load distribution"""
        plt.figure(figsize=(10, 4))
        plt.plot(sensor_x, load / 1e6, 'o-', linewidth=2, markersize=4)
        plt.xlabel('x (m)')
        plt.ylabel('Load (MPa)')
        plt.title('Load Distribution on Top Surface')
        plt.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        
        plt.show()
    
    def plot_training_history(self, history: dict, save_path: Optional[str] = None):
        """Plot training history"""
        plt.figure(figsize=(10, 5))
        plt.semilogy(history['train_loss'], label='Training Loss')
        plt.semilogy(history['val_loss'], label='Validation Loss')
        plt.xlabel('Epoch')
        plt.ylabel('Loss (MSE)')
        plt.title('Training History')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        
        plt.show()


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Main execution function"""
    
    # Set random seeds for reproducibility
    np.random.seed(42)
    torch.manual_seed(42)
    
    # Configuration
    config = {
        'L': 4.0,
        'H': 1.0,
        'E': 210e9,
        'nu': 0.3,
        'n_sensors': 50,
        'n_train_samples': 1000,
        'n_val_samples': 200,
        'n_points_per_sample': 100,
        'batch_size': 32,
        'n_epochs': 200,
        'learning_rate': 1e-3,
        'branch_hidden_dims': [100, 100, 100],
        'trunk_hidden_dims': [100, 100, 100],
        'latent_dim': 100,
    }
    
    print("=" * 80)
    print("DeepONet for 2D Cantilever Beam Linear Elasticity")
    print("=" * 80)
    
    # Generate data
    print("\n1. Generating training data...")
    data_generator = CantileverBeamDataGenerator(
        length=config['L'],
        height=config['H'],
        E=config['E'],
        nu=config['nu'],
        n_sensors=config['n_sensors']
    )
    
    train_loads, train_coords, train_disps = data_generator.generate_dataset(
        config['n_train_samples'],
        config['n_points_per_sample']
    )
    
    val_loads, val_coords, val_disps = data_generator.generate_dataset(
        config['n_val_samples'],
        config['n_points_per_sample']
    )
    
    print(f"   Training samples: {len(train_loads)}")
    print(f"   Validation samples: {len(val_loads)}")
    
    # Create datasets and dataloaders
    train_dataset = CantileverDataset(train_loads, train_coords, train_disps)
    val_dataset = CantileverDataset(val_loads, val_coords, val_disps)
    
    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False)
    
    # Create model
    print("\n2. Creating DeepONet model...")
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    model = DeepONet(
        sensor_count=config['n_sensors'],
        branch_hidden_dims=config['branch_hidden_dims'],
        trunk_hidden_dims=config['trunk_hidden_dims'],
        latent_dim=config['latent_dim'],
        n_outputs=2
    )
    
    # Train model
    print("\n3. Training model...")
    trainer = DeepONetTrainer(model, device=device)
    trainer.train(
        train_loader,
        val_loader,
        n_epochs=config['n_epochs'],
        learning_rate=config['learning_rate']
    )
    
    # Save model
    os.makedirs('models', exist_ok=True)
    trainer.save_model('models/deeponet_cantilever.pth')
    
    # Visualize results
    print("\n4. Visualizing results...")
    visualizer = Visualizer(L=config['L'], H=config['H'])
    
    # Plot training history
    visualizer.plot_training_history(trainer.history, save_path='models/training_history.png')
    
    # Test on a sample
    test_load = data_generator.generate_random_load()
    visualizer.plot_load_distribution(test_load, data_generator.sensor_x, save_path='models/test_load.png')
    visualizer.plot_displacement_field(model, test_load, device=device, save_path='models/test_displacement.png')
    
    print("\n" + "=" * 80)
    print("Training completed successfully!")
    print("Models and visualizations saved to 'models/' directory")
    print("=" * 80)


if __name__ == "__main__":
    main()
