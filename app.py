import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
import mistery # added rts49 2024-04-22 2304z
import streamlit as st

st.title('Hertzsprung–Russell Diagrams and How They Work')

# Read in JWST data
f090w = Table.read('MATCHUP.XYMEEE.F090W', format = 'ascii')
f150w = Table.read('MATCHUP.XYMEEE.F150W', format = 'ascii')

# exclude non-detections
keep = ((f090w['msig'] != 9) & (f150w['msig'] != 9)) & (f090w['ybar'] >= 1000) & (f150w['ybar'] >= 1000)

f090w_keep = f090w[keep]
f150w_keep = f150w[keep]

x_90 = f090w_keep['xbar']
y_90 = f090w_keep['ybar']
x_150 = f150w_keep['xbar']
y_150 = f150w_keep['ybar']

# create color index
cm = f090w_keep['mbar'] - f150w_keep['mbar']

# define x-y zones for each CCD chip
x_hor = [430, 9850]
y_hor = [5730, 4125]

x_vert1 = [2120, 2830]
y_vert1 = [3368, 7423]

x_vert2 = [4690, 5590]
y_vert2 = [2875, 7003]

x_vert3 = [7480, 8270]
y_vert3 = [2483, 6560]

x = [x_hor, x_vert1, x_vert2, x_vert3]
y = [y_hor, y_vert1, y_vert2, y_vert3]


m = []
for i in range(len(y)):
    s = (y[i][1] - y[i][0]) / (x[i][1]-x[i][0])
    m.append(s)
    

#A1
A1 = np.where((y_90 <= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
                (x_90 <= (y_90 - y_vert1[0])/m[2] + x_vert1[0]))

#A2
A2 = np.where((y_90 >= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
                (x_90 <= (y_90 - y_vert1[0])/m[2] + x_vert1[0]))

#A3
A3 = np.where((y_90 <= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 >= (y_90 - y_vert1[0])/m[2] + x_vert1[0]) & \
               (x_90 <= (y_90 - y_vert2[0])/m[1] + x_vert2[0]))

#A4
A4 = np.where((y_90 >= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 >= (y_90 - y_vert1[0])/m[2] + x_vert1[0]) & \
               (x_90 <= (y_90 - y_vert2[0])/m[1] + x_vert2[0]))

#B1
B1 = np.where((y_90 >= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 >= (y_90 - y_vert3[0])/m[2] + x_vert3[0]))

#B2
B2 = np.where((y_90 <= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 >= (y_90 - y_vert3[0])/m[2] + x_vert3[0]))

#B3
B3 = np.where((y_90 >= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 <= (y_90 - y_vert3[0])/m[2] + x_vert3[0]) & \
               (x_90 >= (y_90 - y_vert2[0])/m[1] + x_vert2[0]))

#B4
B4 = np.where((y_90 <= -m[0]*(x_hor[0] - x_90) + y_hor[0]) & \
               (x_90 <= (y_90 - y_vert3[0])/m[2] + x_vert3[0]) & \
               (x_90 >= (y_90 - y_vert2[0])/m[1] + x_vert2[0]))


mbar_90 = f090w_keep['mbar']
mbar_150 = f150w_keep['mbar']

cm_a1 = mbar_90[A1] - mbar_150[A1]
cm_a1adj = cm_a1 - .1

cm_a2 = mbar_90[A2] - mbar_150[A2]
cm_a2adj = cm_a2 - .1

cm_a3 = mbar_90[A3] - mbar_150[A3]
cm_a3adj = cm_a3 - .1

cm_a4 = mbar_90[A4] - mbar_150[A4]
cm_a4adj = cm_a4 - .15

cm_b1 = mbar_90[B1] - mbar_150[B1]

cm_b2 = mbar_90[B2] - mbar_150[B2]
cm_b2adj = cm_b2 - .04

cm_b3 = mbar_90[B3] - mbar_150[B3]
cm_b3adj = cm_b3

cm_b4 = mbar_90[B4] - mbar_150[B4]
cm_b4adj = cm_b4 - .04




st.header('What are Hertzsprung-Russell Diagrams?')

st.markdown('''Hertzsprung-Russell diagrams (often shortened to H-R diagrams) are scatter plots comparing the color index or surface
            temperature of stars to their absolute magnitudes or luminosities. Color index in astronomy refers to a numerical representation
            of the color of a star calculated by finding the difference between two pass band values of a star. The most common color index
            is B-V, which is what will be used here. On this scale, a lower number indicates a bluer and hotter star while a higher number
            indicates a redder and cooler star.''')

st.subheader('Now that we know what an H-R diagram is, how are they useful to us?')

st.markdown('''Through plotting stellar populations on H-R diagrams, there are a number of regions and features that emerge. One of these is
            the diagonal line through the center of the diagram called the main sequence. At this stage in their lives stars are fusing
            hydrogen in their cores and it is the region of an H-R diagram that most stars occupy. Once a star burns all the hydrogen in its
            core, it turns off the main sequence and becomes a giant. From there it begins to burn Helium, and eventually progresses to a
            white dwarf.''')

st.markdown('''The location of a star on an H-R diagram depends on factors such as its age, surface temperature, radius, color index, and
            metallicity. As an example, I have plotted the H-R diagram for the globular cluster M92. Below are a few examples of real stars
            so we can see how they show up on the diagram. You can see that some stars will plot in the very solid main sequence track,
            while others plot in the more sparsely populated region where giants fall.''')


star = st.radio("Select a star",
    ['BAT99-98', "AE Andromedae", "possibly even one more"],
    captions = ["45000 K surface temperature, 37.5 Earth Radii, -0.1 B-V Color Index",
                "20000 K surface temperature, 55 Earth Radii, 0.1 B-V Color Index",
                "blah blah"])

if star == 'BAT99-98':

    fig2, ax = plt.subplots(figsize = (5,5))

    # plot the JWST H-R diagram
    plt.scatter(cm_a1adj, mbar_90[A1], color='#003f5c', s=.5, alpha=.5)
    plt.scatter(cm_a2adj, mbar_90[A2], color='#58508d', s=.5, alpha=.5)
    plt.scatter(cm_a3adj, mbar_90[A3], color='#8a508f', s=.5, alpha=.5)
    plt.scatter(cm_a4adj, mbar_90[A4], color='#bc5090', s=.5, alpha=.5)

    plt.scatter(cm_b1, mbar_90[B1], color='#de5a79', s=.5, alpha=.5)
    plt.scatter(cm_b2adj, mbar_90[B2], color='#ff6361', s=.5, alpha=.5)
    plt.scatter(cm_b3adj, mbar_90[B3], color='#ff8531', s=.5, alpha=.5)
    plt.scatter(cm_b4adj, mbar_90[B4], color='#ffa600', s=.5, alpha=.5)
    
    # plot BAT99-98
    plt.scatter(-0.1, -8.11)

    plt.xlim(-1,1)
    plt.xlabel('B-V Color Index')
    plt.ylim(-15,0)
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    st.pyplot(fig2)
    
elif star == 'AE Andromedae':
    fig3, ax = plt.subplots(figsize = (5,5))

    # plot the JWST H-R diagram
    plt.scatter(cm_a1adj, mbar_90[A1], color='#003f5c', s=.5, alpha=.5)
    plt.scatter(cm_a2adj, mbar_90[A2], color='#58508d', s=.5, alpha=.5)
    plt.scatter(cm_a3adj, mbar_90[A3], color='#8a508f', s=.5, alpha=.5)
    plt.scatter(cm_a4adj, mbar_90[A4], color='#bc5090', s=.5, alpha=.5)

    plt.scatter(cm_b1, mbar_90[B1], color='#de5a79', s=.5, alpha=.5)
    plt.scatter(cm_b2adj, mbar_90[B2], color='#ff6361', s=.5, alpha=.5)
    plt.scatter(cm_b3adj, mbar_90[B3], color='#ff8531', s=.5, alpha=.5)
    plt.scatter(cm_b4adj, mbar_90[B4], color='#ffa600', s=.5, alpha=.5)

    # plot AE Andromedae
    plt.scatter(0.1, -10.2)

    plt.xlim(-1,1)
    plt.xlabel('B-V Color Index')
    plt.ylim(-15,0)
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    st.pyplot(fig3)
    
elif star == 'BAT99-98':
    fig4, ax = plt.subplots(figsize = (5,5))
    plt.scatter(-0.1, -8,11)
    st.pyplot(fig4)




@st.cache_data
def isoPlot(age, metallicity):
    
    # query data for desired age and metallicity from MIST
    iso = mistery.get_isochrone(t=age, FeH=metallicity, photometry='JWST')
    
    # take luminosity and surface temperature data
    iso_lum = iso['log_L']
    iso_temp = iso['log_Teff']
    
    # convert surface temperature to B-V color index and luminosity to absolute magnitude
    # so the MIST isochrone data will align with the JWST data from lab 2
    BV = -0.72 + (7090/(10**iso_temp))
    Mag = -4.83 - 2.5*iso_lum
    
    BV_adj = BV - .54
    Mag_adj = Mag - 3.5
    
    # plot JWST data for background and isochrone
    fig, ax = plt.subplots(figsize = (5,5))

    plt.scatter(cm_a1adj, mbar_90[A1], color='#003f5c', s=.5, alpha=.5)
    plt.scatter(cm_a2adj, mbar_90[A2], color='#58508d', s=.5, alpha=.5)
    plt.scatter(cm_a3adj, mbar_90[A3], color='#8a508f', s=.5, alpha=.5)
    plt.scatter(cm_a4adj, mbar_90[A4], color='#bc5090', s=.5, alpha=.5)

    plt.scatter(cm_b1, mbar_90[B1], color='#de5a79', s=.5, alpha=.5)
    plt.scatter(cm_b2adj, mbar_90[B2], color='#ff6361', s=.5, alpha=.5)
    plt.scatter(cm_b3adj, mbar_90[B3], color='#ff8531', s=.5, alpha=.5)
    plt.scatter(cm_b4adj, mbar_90[B4], color='#ffa600', s=.5, alpha=.5)
    
    plt.plot(BV_adj, Mag_adj)
    
    plt.xlim(-1, 1)
    plt.ylim(-15, 0)
    plt.gca().invert_yaxis()
    plt.xlabel('B-V Color Index')
    plt.ylabel('Magnitude')
    plt.title('Hertzsprung–Russell Diagram')
    st.pyplot(fig)
    
    # return lum

age_slide = st.slider('Age (Gigayears)', 1.0, 19.95, step = 0.5)
met_slide = st.slider('Metallicity', -2.00, 0.50, step = 0.05)

isoPlot(age_slide, met_slide)

