# CEST-LDA-MTRasym

**CEST processing using Lorentzian difference analysis (LDA) and magnetization transfer ratio asymmetry (MTRasym) with inverse Z-spectrum analysis**

Author: Jianpan Huang

Email: jianpanhuang@outlook.com

Affiliation: Department of Diagnostic Radiology, The University of Hong Kong, Hong Kong, China

The demo data is a simulation data created using 5-pool Bloch-McConnell equation with amide fraction varied as 0.0009009, 0.0009009*2, 0.0009009*3, 0.0009009*4. Other parameters remained the same. Therefore, we can see a gradient change in amide (3.5 ppm) map, but not in other CEST maps below.

**You can change the data to your own CEST data, which must include CEST images (img), frequency offsets (offs) and ROI (roi).**

The analysis process include the following procedures:

1. Denoising using MLSVD (optional)
3. DeltaB0 map generation and B0 correction
5. Lorentzian difference analysis (LDA) or magnetization transfer ratio asymmetry (MTRasm)

After running the code, you will see the following fitting process and results:

<img width="1000" alt="image" src="https://github.com/JianpanHuang/CEST-MPLF/assets/43700029/60f6d1e7-0839-4a6b-ab21-ba03b1df7386">

<img width="1003" alt="image" src="https://github.com/JianpanHuang/CEST-LDA-MTRasym/assets/43700029/29303653-c137-4300-9642-9d753ee7ca1d">

If you use the code, please consider citing the references:

[1] Pemmasani Prabakaran R S, Park S W, Lai J H C, et al. Deep‐learning‐based super‐resolution for accelerating chemical exchange saturation transfer MRI. NMR in Biomedicine, 2024: e5130.

[2] Yang Z, Shen D, Chan K W Y, et al. Attention-Based MultiOffset Deep Learning Reconstruction of Chemical Exchange Saturation Transfer (AMO-CEST) MRI. IEEE Journal of Biomedical and Health Informatics, 2024.
