from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors
import os
import glob


def add_rows(grid, row_data: list[dict]) -> None:
	grid.call_api_method('applyTransaction', {'add': row_data})


def style_toc(table_of_contents):
	pass


def fill_toc(table_of_contents, docs: dict[dict]) -> None:
	add_rows(table_of_contents, [entry for entry in docs.values()])


def open_doc_entry(index: int, docs: dict[dict], text_area) -> None:
		text_area.clear()
		with text_area:
			ui.markdown(docs[index]['text'])


def read_docfiles() -> dict[dict]:
	path = os.path.abspath(os.path.dirname(__file__))+'/docs/'
	docfiles = glob.glob(path + '/*.html')

	docs = {}
	for file in docfiles:
		index = ''
		title = ''
		content = []

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


def build(view, **kwargs) -> None:
	view.controller.task_controller.clear_task()
	set_colors()
	header()
	back_button('start', view)

	try:
		docs = read_docfiles()
	except (OSError, PermissionError, IndexError, ValueError) as e:
		print(e)
		view.show_error('could not read docfiles')
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

		style_toc(table_of_contents)
		fill_toc(table_of_contents, docs)

		# ----- text area
		text_area = ui.scroll_area().classes('w-0 pl-14 grow h-full')
		open_doc_entry(1, docs, text_area)

	ui.separator().props('dark')
	footer()

