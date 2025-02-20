The Python 2.7 file DepressionalNetworkv2162023 (last updated February 16, 2023) contains the Depression and DepressionalNetwork classes for 
constructing depressional networks from provided depressional attribute data, and "simulating" time-integrated flow through them, as described in Green and Crumpton (2023).

The Python 2.7 file ExampleNetworkRun provides an example of how the Depression and DepressionalNetowork classes can be used to construct a depressional
network from provided data, simulate the routing of precipitation through the network, and calculate some depressional network properties, as listed and discussed in Green and Crumpton (2023).

The file ExampleNetworkDepressionData.txt in the "Inputs" folder contains the depressional morphology data for the example network featured in Figure 3 in Green and Crumpton (2023).

The file ExampleNetworkJoinTable.txt in the "Inputs" folder contains the network connectivity data for the example network featured in Figure 3 in Green and Crumpton (2023).

The data contained in the ExampleNetworkJoinTable is used to construct the sets of upslope and downslope neighboring depressions for each depression in the network.
This process is performed using the ConstructUpstreamSets and ConstructDownstreamSets methods of the DepressionalNetwork class.

To use the ExampleNetworkRun script download the DepressionalNetworkv2142023 and ExampleNetworkRun python files, and the 
ExampleNetworkDepressionData and ExampleNetworkJoinTable text files to your chosen folders. Change the input and output directories in the ExampleNetworkRun script to point to the ExampleNetworkDepressionData and ExampleNetworkJoinTable files.

The folder named "Outputs" contains copies of the files that are created during execution of the ExampleNetworkRun script. 

Figure 3 in the article is an older, incorrectly scaled version that was accidentally included in the final publication. 
The figure produced by the ExampleNetworkRun script is the correct version.
The article is presently being updated to include the correct version of Figure 3.
