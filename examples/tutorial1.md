**WE ARE POLISHING THIS CODE. IT WILL BE HERE IN THE COMING DAYS - APOLOGIES FOR THE INCONVENIENCE!!**

**Initializing a pipenet object**

_pipenet_ is a Class used to define a vascular network that can be later used to simulated IVIM MRI diffusion signals. 

To initalize the _pipenet_ class the following input parameters are mandatory:

    * nodes:     matrix of node positions (in mm), of size 3 x Nnodes (rows: xpos, ypos, zpos in mm);
                     a node is defined as the input (or as the output) of a pipe
    * radii:     matrix of pipe radii (in mm), of size Nnodes x Nnodes (element (i,j) stores the radius
                 of the pipe connecting node i with node j)
    * qin:       input volumetric flow rate in mm3/s
    * idxin:     index of the input node in the node matrix (nodes[:,idxin] provides the x,y,z coordinates
                 of such a input node)
    * idxout:    index of the output node in the node matrix (nodes[:,idxout] provides the x,y,z coordinates
                 of such an input node)

The easiest way to start building this object is by loading a spreadsheet in .csv format with the following entries:

## Network Data
Here we have a histology image of a mouse kidney slice, where a vascular network has been segmented by tracing visible capillaries.

Each straight capillary segment is annotated [Network_stats.csv](./Network_stats.csv.) with a starting coordinate (Xstart,Ystart,Zstart) and an ending (Xend,Yend,Zend) in µm.
![labels_githhubreoi](https://github.com/user-attachments/assets/0364164f-4f12-4cf2-9fae-3c7f37770e81)

### Data Table

| SegmentNumber | RadiiMean_um      | Xstart      | Ystart     | Zstart | Xend        | Yend       | Zend | Node_at_Xstart_Ystart | Node_at_Xend_Yend |
|---------------|-------------------|-------------|------------|--------|-------------|------------|------|-----------------------|-------------------|
| 1             | 3.3177            | 10056.36    | 10111.85   | 0      | 10083.96    | 10073.92   | 0    | 0                     | 1                 |
| 2             | 2.4370833333333333| 10083.96    | 10073.92   | 0      | 10182.86    | 10164.91   | 0    | 1                     | 2                 |
| 3             | 2.2855166666666666| 10182.86    | 10164.91   | 0      | 10088.94    | 10136.81   | 0    | 2                     | 3                 |
| 4             | 3.306783333333333 | 10056.36    | 10111.85   | 0      | 10088.94    | 10136.81   | 0    | 0                     | 3                 |
| 5             | 1.7060666666666666| 10088.94    | 10136.81   | 0      | 10062.82    | 10167.82   | 0    | 3                     | 4                 |
| 6             | 2.5319333333333334| 10062.82    | 10167.82   | 0      | 10093.35    | 10203.21   | 0    | 4                     | 5                 |
| 7             | 3.3291            | 10093.35    | 10203.21   | 0      | 10141.79    | 10180.75   | 0    | 5                     | 6                 |
| 8             | 2.8676999999999997| 10182.86    | 10164.91   | 0      | 10141.79    | 10180.75   | 0    | 2                     | 6                 |
| 9             | 1.6541833333333333| 10062.82    | 10167.82   | 0      | 10141.79    | 10180.75   | 0    | 4                     | 6                 |
| 10            | 2.0759166666666666| 10141.79    | 10180.75   | 0      | 10183.72    | 10200.81   | 0    | 6                     | 7                 |
| 11            | 2.3742666666666667| 10183.72    | 10200.81   | 0      | 10166.24    | 10247.25   | 0    | 7                     | 8                 |
| 12            | 3.064966666666667 | 10093.35    | 10203.21   | 0      | 10166.24    | 10247.25   | 0    | 5                     | 8                 |
| 13            | 2.6477            | 10166.24    | 10247.25   | 0      | 10185.83    | 10291.86   | 0    | 8                     | 9                 |
| 14            | 1.9695            | 10185.83    | 10291.86   | 0      | 10224.19    | 10278.96   | 0    | 9                     | 10                |
| 15            | 7.66195           | 10224.19    | 10278.96   | 0      | 10215.91    | 10251.86   | 0    | 10                    | 11                |
| 16            | 7.133399999999999 | 10215.91    | 10251.86   | 0      | 10199.99    | 10211.34   | 0    | 11                    | 12                |
| 17            | 1.7906166666666667| 10183.72    | 10200.81   | 0      | 10199.99    | 10211.34   | 0    | 7                     | 12                |
| 18            | 1.3943333333333332| 10199.99    | 10211.34   | 0      | 10270.11    | 10199.48   | 0    | 12                    | 13                |
| 19            | 3.5345333333333335| 10270.11    | 10199.48   | 0      | 10287.6     | 10217.21   | 0    | 13                    | 14                |
| 20            | 6.929283333333333 | 10215.91    | 10251.86   | 0      | 10287.6     | 10217.21   | 0    | 11                    | 14                |
| 21            | 4.208966666666666 | 10224.19    | 10278.96   | 0      | 10278.48    | 10278.52   | 0    | 10                    | 15                |
| 22            | 5.100483333333334 | 10278.48    | 10278.52   | 0      | 10287.6     | 10217.21   | 0    | 15                    | 14                |
| 23            | 2.8586833333333335| 10287.6     | 10217.21   | 0      | 10310.79    | 10215.24   | 0    | 14                    | 16                |
| 24            | 4.654883333333333 | 10265.38    | 10165.38   | 0      | 10287.6     | 10217.21   | 0    | 17                    | 14                |
| 25            | 3.3654833333333333| 10182.86    | 10164.91   | 0      | 10210.56    | 10180.14   | 0    | 2                     | 18                |
| 26            | 4.27815           | 10210.56    | 10180.14   | 0      | 10265.38    | 10165.38   | 0    | 18                    | 17                |
| 27            | 2.386616666666667 | 10265.38    | 10165.38   | 0      | 10271.09    | 10139.0    | 0    | 17                    | 19                |
| 28            | 2.9713666666666665| 10271.09    | 10139.0    | 0      | 10249.72    | 10119.15   | 0    | 19                    | 20                |
| 29            | 2.4316            | 10249.72    | 10119.15   | 0      | 10234.06    | 10135.73   | 0    | 20                    | 21                |
| 30            | 3.2838333333333334| 10234.06    | 10135.73   | 0      | 10210.74    | 10131.01   | 0    | 21                    | 22                |
| 31            | 3.657083333333333 | 10182.86    | 10164.91   | 0      | 10210.74    | 10131.01   | 0    | 2                     | 22                |
| 32            | 2.4389            | 10202.13    | 10106.51   | 0      | 10210.74    | 10131.01   | 0    | 23                    | 22                |
| 33            | 1.8349166666666665| 10202.13    | 10106.51   | 0      | 10172.3     | 10101.79   | 0    | 23                    | 24                |
| 34            | 1.6256833333333331| 10172.3     | 10101.79   | 0      | 10163.87    | 10137.31   | 0    | 24                    | 25                |
| 35            | 2.3390666666666666| 10163.87    | 10137.31   | 0      | 10182.86    | 10164.91   | 0    | 25                    | 2                 |
| 36            | 2.1939666666666664| 10172.3     | 10101.79   | 0      | 10134.48    | 10068.24   | 0    | 24                    | 26                |
| 37            | 1.82385           | 10083.96    | 10073.92   | 0      | 10134.48    | 10068.24   | 0    | 1                     | 26                |


Place the .csv file in the same directory as [script01_initnet.py](./script01_initnet.py), and run the latter to initalize the network.  
