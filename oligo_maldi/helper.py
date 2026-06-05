import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from .sample import Sample
from warnings import deprecated

def key_to_enz_code(key:str, enz_rows, enz_cols) -> int:
    """
    Converts a key to an unambiguous enzyme code based on information in the enz_rows and enz_cols dictionaries.
    """
    # break down key into row and column
    row = key[0]
    col = int(key[1:])
    # select possible enzymes based on row
    row_candidates = []
    for key in enz_rows.keys():
        if row in enz_rows[key]:
            row_candidates.append(key)
    # select possible enzymes based on column
    col_candidates = []
    for key in enz_cols.keys():
        if col in enz_cols[key]:
            col_candidates.append(key)
    # identify enzyme based on row and column
    enz_code = [code for code in row_candidates if code in col_candidates][0]
    return enz_code


def key_to_dntp(key:str, dntp_cols, exclude=[]) -> str:
    """
    Converts a key to an unambiguous dNTP assignment based on information in the dntp_cols dictionary.
    """
    dntp = None # default value to return is None
    key_col = int(key[1:])

    for dntp_key in dntp_cols.keys():
        if key in exclude:
            pass
        else:
            if key_col in dntp_cols[dntp_key]:
                dntp = dntp_key
                break
    return dntp











@deprecated("This was used when processing XML files. The new plain text format does not work with this function.")
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

@deprecated("This was used when processing XML files. The new plain text format does not work with this function.")
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

@deprecated("This was used when working with timsTOF files which were much larger. Unnecessary now.")
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
