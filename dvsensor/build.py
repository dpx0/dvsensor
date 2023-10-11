"""
Run this file inside the same directory as the 'src' directory containing
the DVSensor source code, by running this command in a terminal / command line: python build.py
Make sure you have a working installation of the latest PyInstaller package.
A new 'dist' directory will be created, containing either a 'dvsensor.exe' on Windows or
simply a 'dvsensor' executable binary on Linux. These executables can be run on machines
without having to install python or other dependencies.
"""

import os
import subprocess
import platform
from pathlib import Path
import nicegui

APPNAME = 'dvsensor'
PLATFORM = platform.system()

build_cmd = [
	'python',
	'-m', 'PyInstaller',
	Path(__file__).parent.parent / 'src' / 'main.py',
	'--name', APPNAME,
	'--paths', f'{Path(__file__).parent.parent / "src"}',
	'--hidden-import', 'pages.page_docs',
	'--hidden-import', 'pages.page_input',
	'--hidden-import', 'pages.page_metainf',
	'--hidden-import', 'pages.page_options',
	'--hidden-import', 'pages.page_results',
	'--hidden-import', 'pages.page_start',
	'--add-data', f'{Path(nicegui.__file__).parent}{os.pathsep}nicegui',
	'--add-data', f'{(Path(__file__).parent.parent / "src" / "assets")}{os.pathsep}assets',
	'--add-data', f'{(Path(__file__).parent.parent / "src" / "router_frame.js")}{os.pathsep}.',
	'--add-data', f'{(Path(__file__).parent.parent / "src" / "pages" / "docs")}{os.pathsep}pages/docs'
]

subprocess.call(build_cmd)
