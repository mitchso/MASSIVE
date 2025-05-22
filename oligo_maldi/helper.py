import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from .sample import Sample

def bruker_xml_to_spectra(xml: str):
    """
    Helper function to convert XML files to 'spectra' Elements.
    This works for data pulled from Bruker TIMS-TOF MALDI2 instrument via the export function
    in the DataAnalysis software.
    """
    tree = ET.parse(xml)
    root = tree.getroot()
    ms_spectra = root[2][1][1:]
    return ms_spectra

def spectra_to_samples(spectra: list):
    """
    Helper function to convert 'spectra' from Bruker XMLs to Sample objects.
    """

    samples = []
    for s in spectra:
        well = s.attrib['info'].replace('+MS, ', '')
        peaks = s[0]
        mz_list = [float(peak.attrib['mz']) for peak in peaks]
        intensity_list = [float(peak.attrib['i']) for peak in peaks]

        samples.append(
            Sample(well=well,
                   mz=mz_list,
                   i=intensity_list,
                   noise_cutoff=0.95)
        )

    return samples

def figsave_fast(fig, filename:str):
    """This function will save an image ~4x faster than simply calling plt.savefig
    https://www.scaler.com/topics/matplotlib/save-a-plot-in-matplotlib/
    https://stackoverflow.com/questions/64789437/what-is-the-difference-between-figure-show-figure-canvas-draw-and-figure-canva
    """
    fig.canvas.draw_idle()  # Renders the image without displaying it
    x = np.array(fig.canvas.renderer.buffer_rgba()) # converts the image to bits
    if filename.endswith('.png'):
        return plt.imsave(filename, x)
    else:
        raise ValueError("filename must end with .png")
