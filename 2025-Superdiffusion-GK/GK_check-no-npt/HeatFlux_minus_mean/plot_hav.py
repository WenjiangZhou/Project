import numpy as np
import matplotlib.pyplot as plt
import os


def load_thermal_data(base_dir="."):
    """加载热导率数据"""
    file_path = os.path.join(base_dir, "hac.out")  # 修正文件路径
    
    try:
        # 读取数据，假设至少三列
        data = np.loadtxt(file_path)
        print(f"成功加载数据，形状: {data.shape}")
        
        if data.ndim < 2 or data.shape[1] < 11:  # 需要至少11列，因为要访问索引10
            raise ValueError(f"数据格式不正确: {file_path}，需要至少11列数据")
            
        time = data[:, 0]      # 第一列: 时间
        value = data[:, 10]    # 第11列: kz值
        return time, value
            
    except Exception as e:
        print(f"无法加载数据 {file_path}: {e}")
        return None, None

def plot_thermal_conductivity(base_dir="."):
    """绘制热导率曲线"""
    fig, ax = plt.subplots(figsize=(12, 8))  # 修正括号
    
    time, value = load_thermal_data(base_dir)
    
    if time is not None and value is not None:
        # 绘制热导率曲线
        ax.plot(time, value, 'b-', linewidth=2, label='kz')  # 指定颜色和线宽
        
        # 设置图表属性
        ax.set_title('GPUMD')
        ax.set_xlabel('time (ps)')
        ax.set_ylabel('kz (W/mK)')
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        plt.ylim([-0.5, 3])  # 设置y轴范围
        
        # 保存图表
        #output_file = os.path.join(base_dir, 'thermal_conductivity.png')
        #plt.savefig(output_file, dpi=300, bbox_inches='tight')
        #print(f"图表已保存至: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    # 默认从当前目录读取数据
    plot_thermal_conductivity()
