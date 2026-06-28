"""
Helper functions for MASSIVE.

Warning:
    This section is under active development and is not yet documented. Come back later :)
"""


from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from matplotlib.collections import LineCollection



colour_dict = {'dark grey': '#333333',
               'grey': '#737373',
               'light grey': '#b3b3b3',
               'lighter grey': '#e6e6e6',

               'dark blue': '#007aed',
               'blue': '#5fa8ed',
               'light blue': '#bed6ed',

               'dark red': '#ff2700',
               'red': '#ff7d66',
               'light red': '#ffd4cc',

               'dark green': '#44cc00',
               'green': '#7acc52',
               'light green': '#b1cca3'}

def _savitzky_golay(i: list, window=15, polyorder=3):
    return savgol_filter(i, window_length=window, polyorder=polyorder)


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

def key_to_dntp(key:str, dntp_cols, exclude=()) -> str:
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

def _call_peaks(i_list: list) -> list[int]:
    """
    Uses scipy.signal 'find_peaks' to attempt to call peaks on the graph.
    Returns index positions of the peaks based on the list of intensity values provided.
    """
    peak_indices, properties = find_peaks(i_list,
                                          prominence=200,
                                          height=max(i_list) / 10,
                                          threshold=max(i_list) / 10,
                                          distance=100)

    return list(peak_indices)


def _build_colour_array(x_vals, analyte_ranges, custom_colour_ranges, base_colour, analyte_colour='#d1495b'):
    """
    Returns a per-point colour list the same length as x_vals.

    Priority (lowest → highest):
        base_colour  →  analyte colour  →  custom colours
    """

    # base colour is first applied to everything
    colour_array = [base_colour] * len(x_vals)

    # for each analyte, apply the default analyte colour
    for start, end in analyte_ranges:
        for index, value in enumerate(x_vals):
            if start <= value <= end:
                colour_array[index] = analyte_colour

    # for each custom colour range, apply the custom colour
    for (start, end), custom_colour in custom_colour_ranges.items():
        for index, value in enumerate(x_vals):
            if start <= value <= end:
                colour_array[index] = custom_colour

    return colour_array


def _find_colour_segments(colour_array) -> dict:
    """Takes a colour array and splits it into contiguous segments of the same colour.
    The keys of the dictionary are the segment boundaries, and the values are the colour."""
    colour_segments = {}

    current_colour = colour_array[0]
    segment_start = 0
    for index, colour in enumerate(colour_array):
        if colour != current_colour:
            colour_segments[(segment_start, index+1)] = current_colour
            current_colour = colour
            segment_start = index
        else:
            continue # just do nothing

    # Close the final segment
    colour_segments[(segment_start, len(colour_array))] = current_colour

    return colour_segments


def _plot_segmented(axis, mz, i, colour_array, linewidth=1.5):
    """
    Draws a spectrum as contiguous same-colour runs, emitting one
    ax.plot() call per run. Adjacent runs share a boundary point so
    there are no gaps and no pixel is painted twice.
    """
    colour_segments = _find_colour_segments(colour_array)
    for (seg_start, seg_end), colour in colour_segments.items():
        axis.plot(mz[seg_start:seg_end], i[seg_start:seg_end], c=colour, linewidth=linewidth, label='Experimental')


def _get_nested_attr(object, attr):
    """Helper function for show_sample_positions. Takes a nested attribute string and finds the value.
    i.e. 'enzyme.id' will return the ID of a particular enzyme."""

    attr_list = attr.split(sep='.')

    try:
        # recursively travel through nested objects
        for i, attr in enumerate(attr_list):
            attr_value = getattr(object, attr)
            object = attr_value

        return attr_value

    except AttributeError:
        print(f"Unable to access attribute {attr}. Skipping.")
        return None

def _slice_spectrum(start: int, end: int, mz: list, i: list) -> tuple:
    """
    Slices the spectrum at indicated start and end values, returns a tuple of lists
    """
    mz_slice = []
    i_slice = []
    for j in range(len(mz)):
        if start <= mz[j] <= end:
            mz_slice.append(mz[j])
            i_slice.append(i[j])
    return mz_slice, i_slice