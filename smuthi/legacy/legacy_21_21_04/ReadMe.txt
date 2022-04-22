
Parker Wray edited the direct coupling routine originally written by Amos 
in order to implement faster solutions for calculating direct coupling. 

The direct_coupling.py file in this folder is the original implementation 
written by Amos. 

Parker Wray also added a function to the layers.py to detect a degenerate layer
system. The linear_system.py file then calls this funcitons and checks if degenerate. 
If true, the layer system couplng matrix will be skipped. The layers.py and linear_system.py
files in this folder are the older version written by Amos. 







