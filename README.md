# Basics of Neurons and Network overview
#### A collection of jupyter notebooks built by the Neural Engineering Laboratory at the University of Missouri(Mizzou)

## Running the tutorials

### Google Colab
No local install is needed. Open a folder, open the `.ipynb` file, and click **Open in Colab**. The notebook installs the same pinned packages listed in `requirements.txt` when it detects Colab.

### Local (Poetry)
Requires [Poetry](https://python-poetry.org/docs/#installation), **Python 3.10–3.13**, and a C/C++ toolchain so `nrnivmodl` can compile `.mod` files (NEURON 9.x).

```bash
poetry install
poetry run jupyter notebook
```

`poetry install` is the only local setup command: it installs NEURON, plotting/widget libraries, and Jupyter. Notebooks set `NEURON_MODULE_OPTIONS=-nogui` so NEURON does not abort Jupyter when no GUI display is available. To recompile mechanisms in every tutorial folder:

```bash
poetry run python scripts/compile_all_mods.py
```
### [S1-PassiveMembrane](/S1_PassiveMembrane/)
* #### A simple notebook going over NEURON basics and contains a Hodgkin–Huxley model with widgets
### [S2-ActionPotenial](/S2_ActionPotential/)
* #### An educational notebook going over the fundamentals of how a neuron fires and the voltage gated channels responible for the firing. Contains a NEURON model of a soma with widgets.
### [S3-Burster](/S3_Burster/)
* #### An educational notebook discussing one way neurons can display a bursting effect. Contains a NEURON model of a soma and axon with widgets.
### [S4-Synapses](/S4_Synapses/)
* #### An educational notebook with NEURON model about how synapses operate.
### [S5-CPG](/S5_CPG/)
* #### An educational notebook with a NEURON model about central pattern generates and how differnt movement patterns in horese may arise.
### [S6-STM-WTA Networks](/S6_STM_WTA/)
* #### Educational notebooks talking about short term memory(STM) and winner take all(WTA) networks
### [B1-CreatingSingleCell](/B1_BasicsSingleNeuron/B1_SingleNeuronBio&Model.ipynb)
* #### Goes over the components of a neuron(in particular the soma) and asks the students to make the biological model and electrical circuit. It then goes over how the soma is simulated in NEURON and the default values for properties of the section created. It then goes over inserting the passive and active channels and how to voltage/current clamp for a simulation. Lastly, an interactive simulation is given to show exactly how the soma looks with the bilayer, the fact that it is polarized, and that channels are actually inserted into the soma.
### [B2-CableTheory](/B2_NeuronProperties/B2_Biophysics&Properties.ipynb)
* #### Goes over how to create a dendrite and what electrotonic distance is and how length constant is a part of it. Then it goes over how it is still transient by showing the voltage change in specific segments through the simulation. It then goes over how to calculate the length constant and why it is important.
### [B3-Propagation](/B3_MoreProperties/AdditionalProperties.ipynb)
* #### Goes over what input resistance, the fact that it is a constant, and how to calculate it and asks some questions. It then goes over what the FI curve is and plots it for the given cell from the currents -1nA to 2nA. Lastly, it goes over what Nodes of Ranvier and Myelin sheaths are and how they contribute to AP propagation down an axon by adding them on to the axon. It also asks how certain factors affect the propagation speed.
