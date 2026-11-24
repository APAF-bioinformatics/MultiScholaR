# KIT-183 simulation: outcome

Outcome: **O1**

E1 (camera, 30 checks per cell, upper limit <= 0.075): n = 12 holds (largest upper 0.0737), n = 20 holds (largest upper 0.0711)
E2 at n = 12: gsea holds (largest upper 0.0599); camera_pr holds (largest upper 0.0637)
E3 at n = 12: camera 33877 of 33877 called in the expected direction; gsea 34891 of 34891
Operating characteristic (measured SEs), P(cell passes E1) if every check's true rate were 0.050 / 0.055: n = 12 0.9981 / 0.9587; n = 20 0.9898 / 0.8788

Placeholders:
- R_cam = 6.6 % (C2  proteomics F_ind)
- I_cam = 5.8 to 7.4 % (C2  proteomics F_ind)
- S_cam = 12.8 % (C4b 0.6 proteomics F_none)
- R_gsea_0.3 = 43.0 % (transcriptomics F_main)
- I_gsea_0.3 = 41.9 to 44.1 % (transcriptomics F_main)
- R_cpr_0.3 = 29.9 % (transcriptomics F_main)
- R_cam_6 = 5.0 % (C1 0.1 proteomics F_ind)

Completion:
- n = 6: 800 of 800 completed, 0 did not complete, 0 not started
- n = 12: 800 of 800 completed, 0 did not complete, 0 not started
- n = 20: 800 of 800 completed, 0 did not complete, 0 not started

Determinism: identical
