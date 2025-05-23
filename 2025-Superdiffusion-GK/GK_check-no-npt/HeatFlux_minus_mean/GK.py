import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz
from ase.io import read
import matplotlib.pyplot as plt
import os
from concurrent.futures import ThreadPoolExecutor

# 读取热流数据
def read_heat_flux_data(file_path):
    data = np.loadtxt(file_path, comments='#')
    return data

# 读取初始结构文件并计算体积
def calculate_volume(xyz_file):
    atoms = read(xyz_file)
    return atoms.get_volume()

# 计算自相关函数
def autocorrelation(data, dt, max_lag):
    n = len(data)
    data -= np.mean(data)
    result = np.correlate(data, data, mode='full')[n-1:] / (n - np.arange(n))
    return np.arange(len(result[:max_lag])) * dt, result[:max_lag]

# 计算互相关函数
def cross_correlation(data1, data2, dt, max_lag):
    n = len(data1)
    data1 -= np.mean(data1)
    data2 -= np.mean(data2)
    result = np.correlate(data1, data2, mode='full')[n-1:] / (n - np.arange(n))
    return np.arange(len(result[:max_lag])) * dt, result[:max_lag]

# 计算Green-Kubo热导率和累积热导率
def calculate_kappa(J, v, V, T, dt, max_lag, direction):
    k_B, A_to_m = 1.380649e-23, 1e-10
    V = V * (A_to_m**3)
    #
    J_time, J_corr = autocorrelation(J[:, direction], dt, max_lag)
    v_time, v_corr = autocorrelation(v[:, direction], dt, max_lag)

    cross_time, cross_corr = cross_correlation(J[:, direction], v[:, direction], dt, max_lag)

    #
    kappa_ee = (1 / (k_B * V * T**2)) * cumtrapz(J_corr, J_time, initial=0)
    kappa_ev = (1 / (k_B * V * T**2)) * cumtrapz(cross_corr, J_time, initial=0)
    kappa_vv = (1 / (k_B * V * T**2)) * cumtrapz(v_corr, v_time, initial=0)
    #
    return J_time, J_corr, v_corr, cross_corr, kappa_ee, kappa_ev, kappa_vv 

# 使用最后1/4数据计算数值平均值和标准误差
def calculate_avg_and_error(data):
    quarter_length = len(data[0]) // 2
    avg_values = [np.mean(d[-quarter_length:]) for d in data]
    avg = np.mean(avg_values)
    error = np.std(avg_values) / np.sqrt(len(avg_values))
    return avg, error, avg_values
# 绘图函数
def plot_correlation_and_save(time, data_list, mean_data, title, filename):
    plt.figure(figsize=(4*0.9, 3*0.9))
    #for data in data_list:
    #    plt.plot(time * 1e12, data, color='grey', alpha=0.5)
    plt.semilogy(time * 1e12, mean_data, color='black', linewidth=2, label='Mean')
    plt.xlabel('Time (ps)')
    plt.ylabel('Correlation')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
# 绘图函数
def plot_kee(time, data_list, mean_data, std_data, title, filename):
    plt.figure(figsize=(4*0.9, 3*0.9))
    for data in data_list:
        plt.plot(time * 1e12, data, color='grey', alpha=0.5)
    plt.plot(time * 1e12, mean_data, color='black', linewidth=2, label='Mean')
    plt.plot(time * 1e12, mean_data+std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.plot(time * 1e12, mean_data-std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$\kappa$ (W/mK)')
    plt.ylim([0.2,0.8])
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# 绘图函数
def plot_kev(time, data_list, mean_data, std_data, title, filename):
    plt.figure(figsize=(4*0.9, 3*0.9))
    for data in data_list:
        plt.plot(time * 1e12, data, color='grey', alpha=0.5)
    plt.plot(time * 1e12, mean_data, color='black', linewidth=2, label='Mean')
    plt.plot(time * 1e12, mean_data+std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.plot(time * 1e12, mean_data-std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$\kappa$')
    plt.ylim([1e-8,3e-8])
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# 绘图函数
def plot_kvv(time, data_list, mean_data, std_data, title, filename):
    plt.figure(figsize=(4*0.9, 3*0.9))
    for data in data_list:
        plt.plot(time * 1e12, data, color='grey', alpha=0.5)
    plt.plot(time * 1e12, mean_data, color='black', linewidth=2, label='Mean')
    plt.plot(time * 1e12, mean_data+std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.plot(time * 1e12, mean_data-std_data, color='red', linewidth=1, label='Error', linestyle='--')
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$\kappa$')
    plt.ylim([1e-15,3e-15])
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()



def calculate_y(a, b, c):
    """计算函数 y = a - b²/c 的值（支持数组输入）"""
    return a - b**2 / c

def calculate_y_error(a, b, c, da, db, dc):
    """
    计算函数 y = a - b²/c 的误差（支持数组输入）

    参数:
        a, b, c: 变量数组
        da, db, dc: 变量的误差数组

    返回:
        dy: y的误差数组
    """
    # 计算偏导数（逐元素操作）
    dyda = np.ones_like(a)  # 偏导数 ∂y/∂a = 1
    dydb = -2 * b / c       # 偏导数 ∂y/∂b = -2b/c
    dydc = b**2 / c**2      # 偏导数 ∂y/∂c = b²/c²

    # 计算误差传递（逐元素操作）
    dy = np.sqrt((da * dyda)**2 + (db * dydb)**2 + (dc * dydc)**2)

    return dy


# 主函数
def main(base_dir, xyz_file, T, dt, direction):
    JJ_corr_all, vv_corr_all, Jv_corr_all = [], [], []
    kappa_ee_all, kappa_vv_all, kappa_ev_all = [], [], []
    volume = calculate_volume(xyz_file)

    chunk_size = 200000
    j = read_heat_flux_data('j.txt')
    v1= read_heat_flux_data('v1.txt')
    eV_to_J, AMU_to_kg = 1.60218e-19, 1.66054e-27
    j = j * (eV_to_J**1.5) / (AMU_to_kg**0.5)
    v1 = v1 * (eV_to_J**0.5) * (AMU_to_kg**0.5)
    
    total_data=len(j)

    # 关联时: ps
    t_correlation= 10
    max_lag = round((t_correlation*1e-12)/dt)
    print(max_lag)
    # step, each step=5fs
    #
    for start_idx in range(0, total_data, chunk_size):
        end_idx = min(start_idx + chunk_size, total_data)
        print(f"Processing data chunk: {start_idx} to {end_idx}")
        j_chunk = j[start_idx : end_idx]
        v1_chunk = v1[start_idx : end_idx]

        [J_time, JJ_corr, vv_corr, Jv_corr, kappa_ee, kappa_ev, kappa_vv] = calculate_kappa(j_chunk, v1_chunk, volume, T, dt, max_lag, direction)
        JJ_corr_all.append(JJ_corr)
        vv_corr_all.append(vv_corr)
        Jv_corr_all.append(Jv_corr)
        kappa_ee_all.append(kappa_ee)
        kappa_vv_all.append(kappa_vv)
        kappa_ev_all.append(kappa_ev)


    JJ_corr_mean = np.mean(JJ_corr_all, axis=0)
    vv_corr_mean = np.mean(vv_corr_all, axis=0)
    Jv_corr_mean = np.mean(Jv_corr_all, axis=0)
    kappa_ee_mean = np.mean(kappa_ee_all, axis=0)
    kappa_vv_mean = np.mean(kappa_vv_all, axis=0)
    kappa_ev_mean = np.mean(kappa_ev_all, axis=0)
    kappa_ee_error = np.std(kappa_ee_all, axis=0) / np.sqrt(len(kappa_ee_all))
    kappa_vv_error = np.std(kappa_vv_all, axis=0) / np.sqrt(len(kappa_vv_all))
    kappa_ev_error = np.std(kappa_ev_all, axis=0) / np.sqrt(len(kappa_ev_all))


    kappa_ee_avg, kappa_ee_avg_error, kappa_ee_avg_values = calculate_avg_and_error(kappa_ee_all)
    kappa_vv_avg, kappa_vv_avg_error, kappa_vv_avg_values = calculate_avg_and_error(kappa_vv_all)
    kappa_ev_avg, kappa_ev_avg_error, kappa_ev_avg_values = calculate_avg_and_error(kappa_ev_all)


    plot_correlation_and_save(J_time, JJ_corr_all, JJ_corr_mean, 'JJ Autocorrelation', 'kappa_JJ_corr.pdf')
    plot_correlation_and_save(J_time, vv_corr_all, vv_corr_mean, 'vv Autocorrelation', 'kappa_vv_corr.pdf')
    plot_correlation_and_save(J_time, Jv_corr_all, Jv_corr_mean, 'Jv Autocorrelation', 'kappa_Jv_corr.pdf')
    plot_kee(J_time, kappa_ee_all, kappa_ee_mean, kappa_ee_error, r'$\kappa_{\rm ee}$', 'kappa_ee.pdf')
    plot_kvv(J_time, kappa_vv_all, kappa_vv_mean, kappa_vv_error, r'$\kappa_{\rm vv}$', 'kappa_vv.pdf')
    plot_kev(J_time, kappa_ev_all, kappa_ev_mean, kappa_ev_error, r'$\kappa_{\rm ev}$', 'kappa_ev.pdf')
    
    #print(f"数值平均热导率 (kappa_total_avg ± error): {kappa_total_avg:.6e} ± {kappa_total_avg_error:.6e} W/mK")
    print(f"数值平均基于JJ的热导率 (kappa_ee_avg ± error): {kappa_ee_avg:.6e} ± {kappa_ee_avg_error:.6e} W/mK")
    print(f"数值平均基于vv的热导率 (kappa_vv_avg ± error): {kappa_vv_avg:.6e} ± {kappa_vv_avg_error:.6e} W/mK")
    print(f"数值平均ev的热导率 (kappa_ev_avg ± error): {kappa_ev_avg:.6e} ± {kappa_ev_avg_error:.6e} W/mK")
    print(f"Onsager 热导率 kappa_ee-kappa_ev^2/kappa_vv: {kappa_ee_avg-kappa_ev_avg**2/kappa_vv_avg:.6e} W/mK")

    # 
    np.savetxt('kappa_ee.txt', 
           np.column_stack([J_time*1e12, kappa_ee_mean, kappa_ee_mean+kappa_ee_error, kappa_ee_mean-kappa_ee_error]),
           fmt='%.6e', 
           delimiter='\t', 
           header='time(ps) \tkappa_ee_mean\tkappa_ee_upper\tkappa_ee_lower')
    #
    np.savetxt('kappa_vv.txt',
           np.column_stack([J_time*1e12, kappa_vv_mean, kappa_vv_mean+kappa_vv_error, kappa_vv_mean-kappa_vv_error]),
           fmt='%.6e',
           delimiter='\t',
           header='time(ps) \tkappa_vv_mean\tkappa_vv_upper\tkappa_vv_lower')
    #
    np.savetxt('kappa_ev.txt',
           np.column_stack([J_time*1e12, kappa_ev_mean, kappa_ev_mean+kappa_ev_error, kappa_ev_mean-kappa_ev_error]),
           fmt='%.6e',
           delimiter='\t',
           header='time(ps) \tkappa_ev_mean\tkappa_ev_upper\tkappa_ev_lower')
    # 计算不同关联时间下的热导率
    index_begin = 1
    kappa_total = calculate_y(kappa_ee_mean[index_begin:], kappa_ev_mean[index_begin:], kappa_vv_mean[index_begin:])
    kappa_error = calculate_y_error(kappa_ee_mean[index_begin:], kappa_ev_mean[index_begin:], kappa_vv_mean[index_begin:], kappa_ee_error[index_begin:],kappa_ev_error[index_begin:],kappa_vv_error[index_begin:])
    plt.figure(figsize=(6*0.9, 3*0.9))
    time = J_time[index_begin:]
    plt.plot(time * 1e12, kappa_total, color='black', linewidth=2, label='Mean')
    plt.plot(time * 1e12, kappa_total+kappa_error, color='red', linewidth=1, label='Error', linestyle='--')
    plt.plot(time * 1e12, kappa_total-kappa_error, color='red', linewidth=1, label='Error', linestyle='--')
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$\kappa$ (W/mK)')
    plt.ylim(0.2, 0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig('kappa_total')
    plt.close()
    np.savetxt('kappa_total.txt',
           np.column_stack([time*1e12, kappa_total, kappa_error]),
           fmt='%.6e',
           delimiter='\t',
           header='time(ps) \tkappa_error')



if __name__ == "__main__":
    # direction = x y z
    direction = 2
    main(".", "model.xyz", 650, 1e-15, direction)
