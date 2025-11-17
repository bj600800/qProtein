"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2025/1/13

# Description: Run code for landscape visualization and weblogo generation.
# ------------------------------------------------------------------------------
"""
import configparser
import os
import tempfile
from typing import Dict, Optional

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from weblogo import LogoData, LogoOptions, LogoFormat, eps_formatter, read_seq_data, ColorScheme
from weblogo.color import Color

from qpacking.evolutionary_analysis import function_landscape
from qpacking.common import logger

logger = logger.setup_log(name=__name__)


def plot_landscape(df, input_pos_list, template_name, svg_path):
    fig, ax = plt.subplots(figsize=(18, 14))
    df = df[input_pos_list]
    df = df.dropna(how='all')

    # plotting parameters
    cell_width = 1.2
    cell_height = 1
    padding = 0.5
    edge_padding = 0.3

    def scale_value(value, old_min, old_max, new_min, new_max):
        return new_min + ((value - old_min) / (old_max - old_min)) * (new_max - new_min)

    def get_color(pos_prop):
        if pos_prop > 0.5:
            alpha = scale_value(pos_prop, 0.5, 1, 0, 1)
            color = (1, 0, 0, alpha)

        elif pos_prop < 0.5:
            neg_prop = 1-pos_prop
            alpha = scale_value(neg_prop, 0.5, 1, 0, 1)
            color = (0, 0, 1, alpha)

        else:
            color = (0.5, 0.5, 0.5, 0.3)
        return color

    # cell unit alignment
    for i, amino_acid in enumerate(df.index):
        for j, pos in enumerate(df.columns):
            values = df.at[amino_acid, pos]
            if values is None:
                continue
            pos_prop = values[0]
            color = get_color(pos_prop)
            rect = patches.FancyBboxPatch(
                (j * (cell_width + padding), i * (cell_height + padding)),
                cell_width,
                cell_height,
                boxstyle="round,pad=0.1",
                facecolor=color
            )
            ax.add_patch(rect)

    # plot grid lines
    for i in range(len(df.index) + 1):
        y_pos = i * (cell_height + padding) - padding / 2
        x_start = -padding / 2 - edge_padding
        x_end = len(df.columns) * (cell_width + padding) - padding + edge_padding
        ax.hlines(y_pos, x_start, x_end, colors='black', linewidth=1)

    for j in range(len(df.columns) + 1):
        x_pos = j * (cell_width + padding) - padding / 2
        y_start = -padding / 2 - edge_padding
        y_end = len(df.index) * (cell_height + padding) - padding + edge_padding
        ax.vlines(x_pos, y_start, y_end, colors='black', linewidth=1)

    ax.text(
        -padding / 2 - edge_padding,
        len(df.index) * (cell_height + padding) + cell_height / 2,
        template_name,
        ha='right', va='center',
        fontsize=20,
        fontname='Arial',
        fontweight='bold',
        color='black'
    )

    for j, col in enumerate(df.columns):
        ax.text(
            j * (cell_width + padding) + cell_width / 2,
            len(df.index) * (cell_height + padding) + padding / 2 + edge_padding,
            str(col),
            ha='center', va='center',
            fontsize=20,
            fontname='Arial',
            color='black'
        )

    # row and axis labels
    for i, row in enumerate(df.index):
        ax.text(
            -(2 * edge_padding + padding),
            i * (cell_height + padding) + cell_height / 2,
            str(row),
            ha='center', va='center',
            fontsize=20,
            fontname='Arial',
            color='black'
        )


    ax.set_xlim(-edge_padding, len(df.columns) * (cell_width + padding) - padding + edge_padding)
    ax.set_ylim(len(df.index) * (cell_height + padding) - padding + edge_padding, -edge_padding)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.savefig(svg_path, format='svg')


class CustomColorScheme(ColorScheme):
    """
    A custom color scheme based on a dictionary mapping symbols to colors.
    """
    def __init__(self, logodata: LogoData, scheme: Dict[str, str], default_color: str = "black") -> None:
        """
        Initialize the color scheme with a given mapping of symbols to colors.

        :param logodata: LogoData object containing sequences.
        :param scheme: Dictionary where keys are symbols and values are colors.
        :param default_color: A default color to use when no color is defined for a symbol.
        """
        # Use the parent class constructor to initialize any required fields
        super().__init__(logodata)
        self.scheme = scheme
        self.default_color = Color.from_string(default_color)

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Optional[Color]:
        """
        Return the color for a symbol based on the custom scheme.
        If no color is found, return the default color.
        """
        # If the symbol exists in the custom scheme, return the corresponding color
        if symbol in self.scheme:
            return Color.from_string(self.scheme[symbol])

        return self.default_color


def plot_weblogo(align_seq_df, input_pos_list, eps_path):
    target_res_df = align_seq_df[input_pos_list]
    align_seq_list = []
    for index, row in target_res_df.iterrows():
        align_seq_list.append(''.join(row.tolist()))

    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name
    with open(temp_file_path, 'w') as f:
        f.write('\n'.join(align_seq_list))
    with open(temp_file_path, 'r')as tf:
        seqs = read_seq_data(tf)

    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()

    ####Options start########
    scheme = {
        'K': 'green',
        'R': 'green',
        'H': 'green',
        'D': 'blue',
        'E': 'blue',
        'F': 'purple',
        'W': 'purple',
        'Y': 'purple',
        'M': 'red',
        'V': 'red',
        'A': 'red',
        'L': 'red',
        'I': 'red',
        'G': 'black',
        'N': '#FFB300',
        'Q': '#FFB300',
        'T': '#FFB3FF',
        'S': '#FFB3FF',
        'P': '#FFB3FF',
        'C': '#FFB3FF'

    }
    my_color_scheme = CustomColorScheme(logodata, scheme)

    logooptions.color_scheme = my_color_scheme
    logooptions.show_xaxis = False
    logooptions.show_legend = False
    logooptions.yaxis_label = ''
    logooptions.creator_text = ''
    logooptions.fineprint = ''
    logooptions.text_font = 'Arial-Bold'
    logoformat = LogoFormat(logodata, logooptions)

    ####Options END####

    eps = eps_formatter(logodata, logoformat)
    with open(eps_path, 'wb') as f:
        f.write(eps)
    os.remove(temp_file_path)


def run(label_file, pdb_dir, template_name, input_pos_list, output_xlsx, label_pos_threshold, config_file):
    #### CONFIGURATION PARSER ####
    config = configparser.ConfigParser()
    config.read(config_file)
    usalign_bin = config.get('binary', 'usalign')
    #### END OF CONFIGURATION PARSER ####

    assert label_pos_threshold is float or int, "only int or float is accepted"
    input_pos_list = [int(i) for i in input_pos_list.split(',')]
    figure_dir = os.path.join(pdb_dir, 'figures')
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    svg_path = os.path.join(figure_dir, 'landscape_' + template_name.rstrip('.pdb') + '.svg')
    eps_path = os.path.join(figure_dir, 'weblogo_' + template_name.rstrip('.pdb') + '.eps')

    df_landscape, aligned_seq_df = function_landscape.run(label_file, pdb_dir, template_name, output_xlsx, label_pos_threshold, usalign_bin)
    logger.info("Function landscape results successfully saved to {}".format(output_xlsx))

    plot_landscape(df_landscape, input_pos_list, template_name.rstrip('.pdb'), svg_path)
    logger.info("Landscape plot in vector format (svg) successfully saved to {}".format(svg_path))

    plot_weblogo(aligned_seq_df, input_pos_list, eps_path)
    logger.info("Weblogo plot in vector format (eps) successfully saved to {}".format(eps_path))
