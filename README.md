# REU_summer23

This repository contains two handouts ([sofar_handout](https://github.com/graceguinan/REU_summer23/tree/main/sofar_handout) and [biclustering_handout](https://github.com/graceguinan/REU_summer23/tree/main/biclustering_handout)) and their corresponding codes.  This code was used in [this paper](https://pubs.rsc.org/en/content/articlelanding/2024/va/d3va00361b)[^1], which was published in the Royal Society of Chemistry's Environmental Science: Advances Journal.  A machine learning method called Sparse Orthogonal Factor Regression (SOFAR)[^2] was implemented in both handouts.

The folder [sofar_handout](https://github.com/graceguinan/REU_summer23/tree/main/sofar_handout) demonstrates how SOFAR was used to create a association network between two datasets, [HRMS_data.csv](https://github.com/graceguinan/REU_summer23/blob/main/sofar_handout/data/HRMS_data.csv) and [FDOM_data.csv](https://github.com/graceguinan/REU_summer23/blob/main/sofar_handout/data/FDOM_data.csv).  

The folder [biclustering_handout](https://github.com/graceguinan/REU_summer23/tree/main/biclustering_handout) demonstrates how SOFAR was used as a biclustering technique[^3] on the dataset [HRMS_data.csv](https://github.com/graceguinan/REU_summer23/blob/main/biclustering_handout/data/HRMS_data.csv).

Each folder is structured:

```plaintext
handout/
├── data/
│   └── all datasets used 
├── fit/
│   └── fit created by SOFAR run 
├── helper_functions/
│   └── all helper functions used 
├── images/
│   └── .png images of each plot created  
```




[^1]: Gonsior, M., et al. (2024). Optical Properties and Molecular Differences in Dissolved Organic Matter at the Bermuda Atlantic and Hawai’i ALOHA Time-Series Stations. Environmental Science: Advances 10.1039/D3VA00361B
[^2]: Uematsu, Y., et al. (2019). ”SOFAR: Large-scale association network learning.” IEEE
transactions on information theory 65(8): 4924-4939.
[^3]: Lee, M., Shen, H., Huang, J.Z. and Marron, J.S., 2010. Biclustering via sparse singular
value decomposition. Biometrics, 66(4), pp.1087-1095.



```plaintext
MyProject/
├── docs/
│   ├── index.md
│   └── setup.md
├── src/
│   ├── main.py
│   ├── module/
│   │   └── helper.py
├── tests/
│   ├── test_main.py
│   └── test_helper.py
├── .gitignore
├── README.md
└── requirements.txt
