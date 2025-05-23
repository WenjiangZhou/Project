from pathlib import Path
import sportran as st
import numpy as np
import pandas as pd
from ase.io import read as ase_read, write as ase_write
import matplotlib.pyplot as plt

c = plt.rcParams['axes.prop_cycle'].by_key()['color']

prod_dir = "./"
flux_dir = "./"

Path(flux_dir).mkdir(exist_ok=True, parents=True)

phases = ["alpha"]
T_range = {"alpha": [650,651]}

compute = {}
volumes = {}

i = 0
for phase in phases:
    Ts = np.arange(*T_range[phase], 100)
    print(Ts)
    compute[phase] = {}
    volumes[phase] = {}
    for T in Ts:
        try:
            compute[phase][T] = pd.read_feather(f"./compute.ft")
            print(f"compute.ft loaded.", flush=True)
            volumes[phase][T] = ase_read(
                f"./restart.xyz"
            ).get_volume()
        except FileNotFoundError:
            print(f"./compute.ft NOT found. Loading...", end=" ", flush=True)
            try:
                compute[phase][T] = pd.read_csv(
                    f"{prod_dir}/compute.out",
                    sep=r"\s+",
                    names=[
                        "jv1[1]",
                        "jv2[1]",
                        "jv3[1]",
                        "jv1[2]",
                        "jv2[2]",
                        "jv3[2]",
                        "jv1[3]",
                        "jv2[3]",
                        "jv3[3]",
                        "jk1[1]",
                        "jk2[1]",
                        "jk3[1]",
                        "jk1[2]",
                        "jk2[2]",
                        "jk3[2]",
                        "jk1[3]",
                        "jk2[3]",
                        "jk3[3]",
                        "v1[1]",
                        "v2[1]",
                        "v3[1]",
                        "v1[2]",
                        "v2[2]",
                        "v3[2]",
                        "v1[3]",
                        "v2[3]",
                        "v3[3]",
                    ],
                )

                for ii in range(1, 4):
                    compute[phase][T][f"j[{ii}]"] = (
                        #compute[phase][T][f"jk1[{ii}]"]+
                         compute[phase][T][f"jv1[{ii}]"]
                        #+ compute[phase][T][f"jk2[{ii}]"]
                        + compute[phase][T][f"jv2[{ii}]"]
                        #+ compute[phase][T][f"jk3[{ii}]"]
                        + compute[phase][T][f"jv3[{ii}]"]
                    )
                Path(f"./").mkdir(exist_ok=True, parents=True)
                compute[phase][T].to_feather(f"./compute.ft")
                print("Done. `feather` file saved.", flush=True)
                restart_frame = ase_read(f"{prod_dir}/restart.xyz")
                volumes[phase][T] = restart_frame.get_volume()
                ase_write(f"./restart.xyz", restart_frame)

            except Exception as e:
                print(f"./compute.ft: {e}", flush=True)
        i += 1

for phase in compute:
    for T in compute[phase]:
        # 提取数据
        j_data = compute[phase][T][[f"j[{i}]" for i in range(1, 4)]].to_numpy()
        v1_data = compute[phase][T][[f"v1[{i}]" for i in range(1, 4)]].to_numpy()
        v2_data = compute[phase][T][[f"v2[{i}]" for i in range(1, 4)]].to_numpy()
        total_data = compute[phase][T].shape[0]

upto = 10000

# 假设数据总量是compute[phase][T]['j']的行数
chunk_size = 20000
#total_data = 40000

np.savetxt(
            f"{flux_dir}/j.txt",
            j_data,
            header="j[1] j[2] j[3]",
            fmt="%10.6f",
            delimiter=" "
        )
np.savetxt(
            f"{flux_dir}/v1.txt",
            v1_data,
            header="v1[1] v1[2] v1[3]",
            fmt="%10.6f",
            delimiter=" "
        )
np.savetxt(
            f"{flux_dir}/v2.txt",
            v2_data,
            header="v2[1] v2[2] v2[3]",
            fmt="%10.6f",
            delimiter=" "
        )







#flux = {}
#        compute_fluxes = True
#        if compute_fluxes:
#            flux[phase][T] = st.HeatCurrent(
#                [
#                    compute[phase][T][[f"j[{i}]" for i in range(1, 4)]].to_numpy()[
#                        :upto
#                    ],
#                    compute[phase][T][[f"v1[{i}]" for i in range(1, 4)]].to_numpy()[
#                        :upto
#                    ],
#                    compute[phase][T][[f"v2[{i}]" for i in range(1, 4)]].to_numpy()[
#                        :upto
#                    ],
#                ],
#                DT_FS=5,
#                UNITS="gpumd",
#                TEMPERATURE=T,
#                VOLUME=volumes[phase][T],
#                PSD_FILTER_W=0.025,
#                FREQ_UNITS="THz",
#            )
            #print(volumes[phase][T])
            #print(T)
#            np.save(f"./heatflux.npy", flux[phase][T])
        # Periodogram with given filtering window width
#        j=flux[phase][T]
        #ax = j.plot_periodogram(PSD_FILTER_W=0.4, kappa_units=True, label=r'$\bar{\mathcal{S}}^0_k$')
#        print(j.Nyquist_f_THz)

#        FSTAR_THZ = 25.0
#        jf, ax = j.resample(fstar_THz=FSTAR_THZ, plot=True, freq_units='thz')
        #plt.xlim([0, 20])
       # #ax[1].set_ylim([2, 12]);

        #ax = jf.plot_periodogram(PSD_FILTER_W=0.1)

        #cepstral_analysis
#        jf.cepstral_analysis()

#        results = jf.cepstral_log
#        print(results)

#for phase in compute:
#    flux[phase] = {}
#    for T in compute[phase]:
        # 提取数据
        #j_data = compute[phase][T][[f"j[{i}]" for i in range(1, 4)]].to_numpy()
        #v1_data = compute[phase][T][[f"v1[{i}]" for i in range(1, 4)]].to_numpy()
        #v2_data = compute[phase][T][[f"v2[{i}]" for i in range(1, 4)]].to_numpy()

        # 保存到文本文件
        #np.savetxt(
        #    f"{flux_dir}/j.txt",
        #    j_data,
        #    header="j[1] j[2] j[3]",
        #    fmt="%10.6f",
        #    delimiter=" "
        #)
        #np.savetxt(
        #    f"{flux_dir}/v1.txt",
        #    v1_data,
        #    header="v1[1] v1[2] v1[3]",
        #    fmt="%10.6f",
        #    delimiter=" "
        #)
        #np.savetxt(
        #    f"{flux_dir}/v2.txt",
        #    v2_data,
        #    header="v2[1] v2[2] v2[3]",
        #    fmt="%10.6f",
        #    delimiter=" "
        #)

        #compute_fluxes = True
        #if compute_fluxes:
        #    flux[phase][T] = st.HeatCurrent(
        #        [
        #            compute[phase][T][[f"j[{i}]" for i in range(1, 4)]].to_numpy()[
        #                :upto
        #            ],
        #            compute[phase][T][[f"v1[{i}]" for i in range(1, 4)]].to_numpy()[
        #                :upto
        #            ],
        #            compute[phase][T][[f"v2[{i}]" for i in range(1, 4)]].to_numpy()[
        #                :upto
        #            ],
        #        ],
        #        DT_FS=5,
        #        UNITS="gpumd",
        #        TEMPERATURE=T,
        #        VOLUME=volumes[phase][T],
        #        PSD_FILTER_W=0.025,
        #        FREQ_UNITS="THz",
        #    )
            #print(volumes[phase][T])
            #print(T)
        #    np.save(f"./heatflux.npy", flux[phase][T])
        # Periodogram with given filtering window width
        #j=flux[phase][T]
        #ax = j.plot_periodogram(PSD_FILTER_W=0.4, kappa_units=True, label=r'$\bar{\mathcal{S}}^0_k$')
        #print(j.Nyquist_f_THz)
        # compare with the spectrum of the energy flux
        #jen = st.HeatCurrent(
        #    [
        #        compute[phase][T][[f"j[{i}]" for i in range(1, 4)]].to_numpy()[
        #            :upto],
        #        ],
        #        DT_FS=5,
        #        UNITS="gpumd",
        #        TEMPERATURE=T,
        #        VOLUME=volumes[phase][T],
        #        PSD_FILTER_W=0.2,
        #        FREQ_UNITS="THz",
        #    )

        #j = jen
        #ax = jen.plot_periodogram(axes=ax, PSD_FILTER_W=0.2, kappa_units=True, label=r'$\mathcal{S}^0_k$')
        #plt.xlim([0, 20])
        #ax[0].set_ylim([0, 1]);
        #ax[1].set_ylim([2, 12]);
        #ax[0].legend(); ax[1].legend();

        #resample
        #FSTAR_THZ = 25.0
        #jf, ax = j.resample(fstar_THz=FSTAR_THZ, plot=True, freq_units='thz')
        #plt.xlim([0, 20])
        #ax[1].set_ylim([2, 12]);

        #ax = jf.plot_periodogram(PSD_FILTER_W=0.1)

        #cepstral_analysis
        #jf.cepstral_analysis()

        # Cepstral Coefficients
        #print('c_k = ', jf.cepf.logpsdK)

        #ax = jf.plot_ck()
        #ax.set_xlim([0, 25])
        #ax.set_ylim([-0.2, 1.0])
        #ax.grid();

        # AIC function
        #f = plt.figure()
        #plt.plot(jf.cepf.aic, '.-', c=c[0])
        #plt.xlim([0, 50])
        #plt.ylim([14600, 14800]);

        #print('K of AIC_min = {:d}'.format(jf.cepf.aic_Kmin))
        #print('AIC_min = {:f}'.format(jf.cepf.aic_min))
        # 清空当前图形窗口的内容
        #plt.close()

        # L_0 as a function of cutoff K
        #ax = jf.plot_L0_Pstar()
        #ax.set_xlim([0, 50])
        #ax.set_ylim([14.75, 16.5]);

        #print('K of AIC_min = {:d}'.format(jf.cepf.aic_Kmin))
        #print('AIC_min = {:f}'.format(jf.cepf.aic_min))

        # kappa as a function of cutoff K
        #ax = jf.plot_kappa_Pstar()
        #ax.set_xlim([0, 50])
        #ax.set_ylim([0, 1.0]);

        #print('K of AIC_min = {:d}'.format(jf.cepf.aic_Kmin))
        #print('AIC_min = {:f}'.format(jf.cepf.aic_min))

        #results = jf.cepstral_log
        #print(results)

        # filtered log-PSD
        #ax = j.plot_periodogram(0.5, kappa_units=True)
        #ax = jf.plot_periodogram(0.5, axes=ax, kappa_units=True)
        #ax = jf.plot_cepstral_spectrum(axes=ax, kappa_units=True)
        #ax[0].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
        #ax[1].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
        #plt.xlim([0, 20])
        #ax[0].set_ylim([0, 0.8]);
        #ax[1].set_ylim([4, 6]);
        #ax[0].legend(['original', 'resampled', 'cepstrum-filtered'])
        #ax[1].legend(['original', 'resampled', 'cepstrum-filtered']);

        #plt.show()


