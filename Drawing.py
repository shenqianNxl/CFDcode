import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def read_and_plot_fluid_data():
    # 读取CSV文件
    df = pd.read_csv('fluid_data.csv')
    
    # 提取数据
    X = df['X'].values
    Y = df['Y'].values
    rho = df['rho'].values
    u = df['u'].values
    v = df['v'].values
    P = df['P'].values
    
    # 创建三角剖分用于等高线绘图
    triang = tri.Triangulation(X, Y)
    
    # 设置图形
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Fluid Dynamics Simulation Results', fontsize=16)
    
    # 绘制密度分布
    contour1 = axes[0,0].tricontourf(triang, rho, levels=15, cmap='viridis')
    axes[0,0].set_title('Density Distribution')
    axes[0,0].set_xlabel('X')
    axes[0,0].set_ylabel('Y')
    axes[0,0].set_aspect('equal')
    axes[0,0].set_xlim(0, 4)
    axes[0,0].set_ylim(0, 2.4)
    plt.colorbar(contour1, ax=axes[0,0])
    
    # 绘制x方向速度分布
    contour2 = axes[0,1].tricontourf(triang, u, levels=15, cmap='plasma')
    axes[0,1].set_title('Velocity in X Direction')
    axes[0,1].set_xlabel('X')
    axes[0,1].set_ylabel('Y')
    axes[0,1].set_aspect('equal')
    axes[0,1].set_xlim(0, 4)
    axes[0,1].set_ylim(0, 2.4)
    plt.colorbar(contour2, ax=axes[0,1])
    
    # 绘制y方向速度分布
    contour3 = axes[1,0].tricontourf(triang, v, levels=15, cmap='plasma')
    axes[1,0].set_title('Velocity in Y Direction')
    axes[1,0].set_xlabel('X')
    axes[1,0].set_ylabel('Y')
    axes[1,0].set_aspect('equal')
    axes[1,0].set_xlim(0, 4)
    axes[1,0].set_ylim(0, 2.4)
    plt.colorbar(contour3, ax=axes[1,0])
    
    # 绘制压力分布
    contour4 = axes[1,1].tricontourf(triang, P, levels=15, cmap='hot')
    axes[1,1].set_title('Pressure Distribution')
    axes[1,1].set_xlabel('X')
    axes[1,1].set_ylabel('Y')
    axes[1,1].set_aspect('equal')
    axes[1,1].set_xlim(0, 4)
    axes[1,1].set_ylim(0, 2.4)
    plt.colorbar(contour4, ax=axes[1,1])
    
    # 调整布局
    plt.tight_layout()
    plt.savefig('fluid_dynamics_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 打印数据统计信息
    print("Data Statistics:")
    print(f"Density: min={rho.min():.4f}, max={rho.max():.4f}, mean={rho.mean():.4f}")
    print(f"Velocity X: min={u.min():.4f}, max={u.max():.4f}, mean={u.mean():.4f}")
    print(f"Velocity Y: min={v.min():.4f}, max={v.max():.4f}, mean={v.mean():.4f}")
    print(f"Pressure: min={P.min():.4f}, max={P.max():.4f}, mean={P.mean():.4f}")

if __name__ == "__main__":
    read_and_plot_fluid_data()