# CalcCEA2034

CalcCEA2034 is a simple tool to calculate CEA2034 Spinorama data as CSV text files from raw horizontal and vertical measurements (e.g. as provided by Audio Science Review loudspeaker reviews). Please note that the raw measurement data may be published under some license which may effect the calculated data as well. E.g. Audio Science Review data is currently published under Creative CommonsAttribution-NonCommercial-ShareAlike 4.0 International License. ChatGPT AI was used to generate CalcCEA2034.py.

The calculated data (CEA2034, Early Refections and Estimated In-Room Response) can directly be imported into various tools, e.g. [REW](https://www.roomeqwizard.com/), [FreqRespGraph](https://github.com/rwhomeaudio/FreqRespGraph) or spreadsheet software for display and further processing.


# Installation

## From Binary Release
For Windows 11 and Ubuntu 24.04 a [PyInstaller](https://pyinstaller.org/en/stable/) standalone onefile binary build is provided in the [Releases](https://github.com/rwhomeaudio/CalcCEA2034/releases). Download it and extract the executable. If it doesn't work on your system or a different platform is required install it from source.
## From Source
1. Download and install Python 3.X on your system: [https://www.python.org/downloads/](https://www.python.org/downloads/)
2. Open a shell / command prompt and install numphy:
```
pip install numpy
```
3. Optionally install git command line client: [https://github.com/git-guides/install-git](https://github.com/git-guides/install-git)
4. Download CalcCEA2034:
   - Clone git repository: `git clone https://github.com/rwhomeaudio/CalcCEA2034` or
   - Download ZIP archive from [https://github.com/rwhomeaudio/CalcCEA2034](https://github.com/rwhomeaudio/CalcCEA2034) using "Code->Download ZIP" button.
5. Test the installation by running: `python CalcCEA2034.py -h`
# Usage
```
python CalcCEA2034.py -h
usage: CalcCEA2034.py [-h] [--output-dir OUTPUT_DIR] [input_dir]

Calculate CEA-2034 data curves from Klippel horizontal and vertical SPL exports.

positional arguments:
  input_dir             Folder containing "SPL Horizontal.txt" and "SPL Vertical.txt" (default: current folder).

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
                        Output folder (default: same as input folder).

Example:

  Download measurement ZIP archive from an ASR review (e.g. https://www.audiosciencereview.com/forum/index.php?threads/ascend-cbm-170se2-speaker-review.72683/),
  extract files and run

  CalcCEA2034.py

  It creates

    CEA2034/Early Reflections DI.txt
    CEA2034/On Axis.txt
    CEA2034/Early Reflections.txt
    CEA2034/Sound Power.txt
    CEA2034/DI offset.txt
    CEA2034/Listening Window.txt
    CEA2034/Sound Power DI.txt
    Estimated In-Room Response.txt
    Early Reflections/Front Wall Bounce.txt
    Early Reflections/All.csv
    Early Reflections/Total Early Reflection.txt
    Early Reflections/Floor Bounce.txt
    Early Reflections/Rear Wall Bounce.txt
    Early Reflections/Ceiling Bounce.txt
    Early Reflections/Side Wall Bounce.txt

  from the raw measurement data ("SPL Horizontal.txt", "SPL Vertical.txt").
```
# Examples
