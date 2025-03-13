import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
params = {'legend.fontsize': 11,
          'legend.handlelength': 0.5, 'xtick.labelsize':11, 'ytick.labelsize':11, 
          'font.family':'arial','legend.labelspacing':0.2, 'legend.title_fontsize':12, 
          'legend.borderpad':0.2, 'legend.frameon':False, 'legend.columnspacing':1,'pdf.fonttype':42}
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update(params)
import pyprocar
from pymatgen.io.vasp.outputs import Outcar

def plot_bands(PROCAR_file,OUTCAR_file,plot_name,energy_range = [-2.5,5], print_energies = True ): 
    '''
    Plots bands/orbitals from VASP output files from a spin polarized DFT calcualtions. Colors orbitals based on d character. 
    Inputs: 
        PROCAR_file = PROCAR file from a spin polarized VASP calcualtion 
        OUTCAR_file = OUTCAR files from a spin polarized VASP calcualtion
        plot_name: files name to save the final plot output at 
        energy_range: energy range for the orbtials you wish to plot. Range is in eV and relative to the energy of the highest occupied molecular obrital
        print_energies: if True energies of the HOMO, band gap, and spin flip gaps are printed
    '''
    #read in and parse the PROCAR file
    parse = pyprocar.ProcarParser()
    parse.readFile(PROCAR_file)  
    orbital_num = int(parse.bands.shape[1]/2) #number of orbitals of each spin polarization
    
    #seperate energies of bands/orbirals by spin up and spin down
    bands_up_iso_0MV = parse.bands[0][0:orbital_num]
    bands_down_iso_0MV = parse.bands[0][orbital_num:] 
    
    #find the character of the bands/orbitals
    spd_iso_0MV = parse.spd 
    fermi = pyprocar.UtilsProcar()
    
    #Get the Fermi energy
    fermi_energy=fermi.FermiOutcar(OUTCAR_file) 
    
    # find the HOMO index and energy via the Fermi energy
    HOMO_indx = sum(bands_up_iso_0MV <fermi_energy)-1
    HOMO_energy_0MV = bands_up_iso_0MV[HOMO_indx] 
    if print_energies == True:
        print('HOMO energy',HOMO_energy_0MV) 
        print('band_gap',bands_up_iso_0MV[HOMO_indx+1]-bands_up_iso_0MV[HOMO_indx])
        print('spin flip gap',bands_down_iso_0MV[HOMO_indx-1]-bands_up_iso_0MV[HOMO_indx])
        print('spin flip gap',bands_down_iso_0MV[HOMO_indx+1]-bands_up_iso_0MV[HOMO_indx])
        print('Final Energy',Outcar(OUTCAR_file).final_energy)
        
    # Find the character of the orbital/bands
    spd_iso_norm_0MV = np.divide(spd_iso_0MV[0,:,0,:,10], spd_iso_0MV[0,:,0,-1,10][:,None])

    #plot the data
    fig, ax = plt.subplots()
    x_1 = np.zeros(len(bands_up_iso_0MV) )
    x_2 = np.zeros(len(bands_down_iso_0MV) )+0.01
    fig.set_figwidth(3)
    fig.set_figheight(5)
    plt.scatter(x_1, bands_up_iso_0MV-HOMO_energy_0MV, s=1000, marker="_", linewidth=3, zorder=3, label='d$_{z^{2}}$',c=spd_iso_norm_0MV[:96,0], cmap='viridis',vmin=0, vmax=1) 
    plt.scatter(x_2, bands_down_iso_0MV-HOMO_energy_0MV, s=1000, marker="_", linewidth=3, zorder=3, label='d$_{z^{2}}$',c=spd_iso_norm_0MV[96:,0], cmap='viridis',vmin=0, vmax=1) 
    ax.axhline(y=0, linestyle='dashed', c='k')
    plt.ylim(energy_range) 
    plt.xlim([-0.005,0.015]) 
    ax.set_xticklabels(['','Spin \n Up','Spin \n Down'],size=11)
    plt.ylabel('Energies (eV)',size=14)
    plt.colorbar() 
    plt.savefig(plot_name,dpi=400)