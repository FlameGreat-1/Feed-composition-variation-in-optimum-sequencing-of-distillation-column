# Feed Composition Variation in Optimum Sequencing of Distillation Column

## Description

This project implements a graphical user interface (GUI) for analyzing and optimizing the sequencing of distillation columns based on feed composition variations. 
It uses PyQt5 for the GUI and incorporates complex calculations for distillation column design and cost estimation.

## Features

- Interactive GUI for input of feed flow rate and composition
- Calculation of optimal distillation column sequences
- Detailed output for each column in the sequence, including:
  - Design parameters (diameter, height, number of plates)
  - Cost breakdown (steam, cooling water, capital, operating costs)
  - Heat exchanger duties and areas
- Summary of total costs for the entire sequence
- Additional design information for each column

## Requirements

- Python 3.6+
- PyQt5
- NumPy

## Installation

1. Clone this repository:
   git clone https://github.com/FlameGreat-1/distillation-column-Algorithm.git


2. Navigate to the project directory:
   cd distillation-column-algorithm

3. Install the required packages:
   pip install -r requirements.txt


## Usage

1. Run the main script:
   python DC_Sequence_Algo.py       ## Without Graphical User Interface
   python Algorithm_GUI.py          ## For Graphical User Interface


2. In the GUI:
- Enter the feed flow rate in the first input field
- Enter the feed composition (comma-separated values) in the second input field
- Click the "Compute" button to run the calculations

3. View the results in the tabbed output area:
- Each "Column" tab shows details for a specific column in the sequence
- The "Cost Summary" tab provides an overview of total costs
- The "Distillation Column Design Information" tab shows additional design parameters

## Contributing

Contributions to this project are welcome. Please fork the repository and submit a pull request with your changes.

## License

[MIT]

## Contact

[Flame Great] - Contact me @:
               -[eugochukwu77@gmail.com]
               - [12348136872013]

Project Link: https://github.com/FlameGreat-1/Feed-composition-variation-in-optimum-sequencing-of-distillation-column

