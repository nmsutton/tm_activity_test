int Axo_Axonic=sim.createGroup("Axo_Axonic", 
                            1, INHIBITORY_NEURON, ANY, GPU_CORES);  
int Pyramidal=sim.createGroup("Pyramidal", 
                            1, EXCITATORY_NEURON, ANY, GPU_CORES);

sim.setNeuronParameters(Axo_Axonic, 165.0f, 0.0f, 3.961462878f, 0.0f, -57.09978287f, 0.0f, -51.71875628f, 
                                0.0f, 0.004638608f, 0.0f, 8.683644937f, 0.0f, 27.79863559f, 0.0f, -73.96850421f, 0.0f, 
                                15.0f, 0.0f, 1); // C,k,vr,vt,a,b,vpeak,c,d
sim.setNeuronParameters(Pyramidal, 366.0f, 0.0f, 0.792338703789581f, 0.0f, -63.2044008171655f, 0.0f, -33.6041733124267f, 
                                0.0f, 0.00838350334098279f, 0.0f, -42.5524776883928f, 0.0f, 35.8614648558726f, 0.0f, -38.8680990294091f, 0.0f, 
                                588.0f, 0.0f, 1);

/* neuron connection parameters */
sim.connect(Axo_Axonic, Pyramidal, "one-to-one", 1.0f, 1.0f, 
    RangeDelay(1), RadiusRF(-1), SYN_FIXED, 3.644261648, 0.0);

/* STP parameters */
sim.setSTP(Axo_Axonic, Pyramidal, true, STPu(0.259914361, 0.0f),
                                     STPtauU(17.20004939, 0.0f),
                                     STPtauX(435.8103009, 0.0f),
                                     STPtdAMPA(10.71107251, 0.0f),
                                     STPtdNMDA(150.0, 0.0f),
                                     STPtdGABAa(10.71107251, 0.0f),
                                     STPtdGABAb(150.0, 0.0f),
                                     STPtrNMDA(0.0f, 0.0f),
                                     STPtrGABAb(0.0f, 0.0f));