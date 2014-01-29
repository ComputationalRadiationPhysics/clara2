'''
  Copyright 2014 Richard Pausch
 
  This file is part of Clara 2.
 
  Clara 2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  Clara 2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with Clara 2.
  If not, see <http://www.gnu.org/licenses/>.
'''




from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


omega_max = 7.0e19 #max frequency   
theta_max = 2.5 #max angle norm. to gamma
tickfontsize=18 #plot font
labelfontsize=28 #plot font

plt.figure(0, figsize=(15, 15))
SP = plt.subplot(111, autoscale_on=False, xlim=(0, omega_max), ylim=(0, theta_max))

phi_index = 0

for id_string in ["phi_0", "phi_90"]:
    plt.title("Radiation Undulator", fontsize=labelfontsize)
    plt.xlabel(r"$\omega$", fontsize=28)
    plt.xticks(size=tickfontsize)
    plt.ylabel(r"$\theta / \gamma$", fontsize=labelfontsize)
    plt.yticks(size=tickfontsize)
    data = loadtxt("my_spectrum_all_%03d.dat" %phi_index)
    plt.imshow(data, extent=(0, omega_max,  0, theta_max), interpolation='nearest', aspect='auto', origin='lower')
    if id_string == "phi_0":
    	CB = plt.colorbar()
        CB.set_label(r"$\frac{\mathrm{d} ^2I}{\mathrm{d} \omega \mathrm{d} \Omega}$", fontsize=labelfontsize)

    #else:
    #    CB.update_bruteforce()
    #plt.show()
    #plt.draw()
    plt.savefig("plot_omega_theta_"+id_string+".pdf", format='pdf')
    SP.clear()
    phi_index = phi_index+1

print "Done with drawing the results"
