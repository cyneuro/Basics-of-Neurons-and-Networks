# Computional Neuroscience Tutorials
#### A collection of jupyter notebooks builts by the Neural Engineering Laboratory at the University of Missouri(Mizzou)

## All tutorials can be on either a local computer or on Google Colab. To run on Colab there is no need to download the repo and instead click on the folder and then the .ipynb file. There will be a button at the start of the notebook which says 'Open in Colab'. Click this and the notebook will loaded into Colab. If running on a local computer machine sure that Jupyter is installed with the python packages neuron and ipywidgets.
### [S1-PassiveMembrane](/S1_PassiveMembrane/)
* #### A simple notebook going over NEURON basics and contains a Hodgkin–Huxley model with widgets
### [S2-ActionPotenial](/S2_ActionPotential/)
* #### An educational notebook going over the fundamentals of how a neuron fires and the voltage gated channels responible for the firing. Contains a NEURON model of a soma with widgets.
### [S3-Burster](/S3_Burster/)
* #### An educational notebook discussing one way neurons can display a bursting effect. Contains a NEURON model of a soma and axon with widgets.
### [S4-Synapses](/S4_Synapses/)
* #### An educational notebook with NEURON model about how synapses operate.
### [S5-CPG](/S5_CPG/)
* #### An educational notebook with a NEURON modle about central pattern generates and how differnt movement patterns in horese may arise.
### [S6-STM-WTA Networks](/S6_STM_WTA/)
* #### Educational notebooks talking about short term memory(STM) and winner take all(WTA) networks
### [B1-CreatingSingleCell](/B1_BasicsSingleNeuron/B1_SingleNeuronBio&Model.ipynb)
* #### Goes over the components of a neuron(in particular the soma) and asks the students to make the biological model and electrical circuit. It then goes over how the soma is simulated in NEURON and the default values for properties of the section created. It then goes over inserting the passive and active channels and how to voltage/current clamp for a simulation. Lastly, an interactive simulation is given to show exactly how the soma looks with the bilayer, the fact that it is polarized, and that channels are actually inserted into the soma.
### [B2-CableTheory](/B2_NeuronProperties/B2_Biophysics&Properties.ipynb)
* #### Goes over how to create a dendrite and what electrotonic distance is and how length constant is a part of it. Then it goes over how it is still transient by showing the voltage change in specific segments through the simulation. It then goes over how to calculate the length constant and why it is important.
### [B3-Propagation](/B3_MoreProperties/AdditionalProperties.ipynb)
* #### Goes over what input resistance, the fact that it is a constant, and how to calculate it and asks some questions. It then goes over what the FI curve is and plots it for the given cell from the currents -1nA to 2nA. Lastly, it goes over what Nodes of Ranvier and Myelin sheaths are and how they contribute to AP propagation down an axon by adding them on to the axon. It also asks how certain factors affect the propagation speed.