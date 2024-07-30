[Network_stats_example.csv](https://github.com/user-attachments/files/16426442/Network_stats_example.csv)[Network_stats_example.csv](https://github.com/user-attachments/files/16426422/Network_stats_example.csv)**WE ARE POLISHING THIS CODE. IT WILL BE HERE IN THE COMING DAYS - APOLOGIES FOR THE INCONVENIENCE!!**

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

[Uploading NetSegmentNumber,RadiiMean_um,Xstart,Ystart,Zstart,Xend,Yend,Zend,Node_at_Xstart_Ystart,Node_at_Xend_Yend
1,4.13,16841.83,27271.38,0,16811.28,27324.12,0,0,1
2,3.71,16841.83,27271.38,0,16889.67,27271.85,0,0,2
work_stats_example.csv…]()

Place the .csv file in the same directory as script01_initnet.py, and run the latter to initalize the network.  
