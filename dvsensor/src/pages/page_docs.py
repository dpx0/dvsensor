import logging
from pathlib import Path
from nicegui import ui
from functools import lru_cache
from typing import Any
from interface import Controller
from .pagebuilder import page_builder
from .ui_elements import back_button, show_error_page


def add_rows(grid, row_data: list[dict]) -> None:
	grid.call_api_method('applyTransaction', {'add': row_data})


def fill_toc(table_of_contents, docs: dict[int, dict]) -> None:
	add_rows(table_of_contents, [entry for entry in docs.values()])


def open_doc_entry(index: int, docs: dict[int, dict], text_area) -> None:
	text_area.clear()
	with text_area:
		if index in docs.keys():
			ui.markdown(docs[index]['text'])


@lru_cache()
def read_docfiles() -> dict[int, dict]:
	docfiles = (Path(__file__).parent / 'docs').glob('*.html')
	docs = {}

	for file in docfiles:
		index: int | None = None
		title: str = ''
		content: list[str] = []

		with open(file) as f:
			while True:
				line = f.readline()
				if not line:
					break
				if line.startswith('<-- Index: '):
					index = int(line.split('<-- Index:')[1].split('-->')[0].strip())
				elif line.startswith('<-- Title: '):
					title = line.split('<-- Title:')[1].split('-->')[0].strip()
				else:
					content.append(line)
		docs[index] = {'index': index, 'title': title, 'text': ''.join(content)}
	return docs


@page_builder()
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	back_button(controller, '/')

	try:
		docs = read_docfiles()
	except (OSError, PermissionError, IndexError, ValueError) as e:
		logging.error(f'could not read docfiles: {e}')
		show_error_page('could not read docfiles', controller)
		return

	ui.separator().props('dark').classes('mt-4')
	with ui.grid(columns=2).classes('w-full flex mt-6').style('height: 74vh'):

		# ----- table of contents
		table_of_contents = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Title', 'field': 'title'},
				{'headerName': 'index', 'field': 'index', 'hide': True}
			],
			'rowData': [],
			'rowSelection': 'multiple',
		}, theme='alpine-dark')\
			.classes('w-1/5 h-full')\
			.on('rowSelected',
				lambda event: open_doc_entry(event.args['data']['index'], docs, text_area)
				if event.args['selected'] else None)

		fill_toc(table_of_contents, docs)

		# ----- text area
		text_area = ui.scroll_area().classes('w-0 pl-14 grow h-full')
		open_doc_entry(1, docs, text_area)

	ui.separator().props('dark')
