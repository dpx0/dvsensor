from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors
import os
import glob


def find_docfiles() -> list[str]:
	path = os.path.abspath(os.path.dirname(__file__))+'/docs/'
	files = glob.glob(path + '/*.md')
	return files


def read_docfiles(docfiles: list[str]) -> list[tuple[int, str, str]]:

	docs = []
	for file in docfiles:
		index = ''
		title = ''
		content = []

		with open(file) as f:
			while True:
				line = f.readline()
				if not line: break
				if line.startswith('<> Index: '):
					index = line.split('<> Index: ')[1].strip()
				elif line.startswith('<> Title: '):
					title = line.split('<> Title: ')[1].strip()
				else:
					content.append(line)
		docs.append((index, title, ''.join(content)))
	return docs


def build(view, **kwargs) -> None:
	view.controller.clear_job()
	set_colors()
	header()
	back_button('start', view)

	docs = read_docfiles(find_docfiles())
	print(docs)

	with ui.grid(columns=2).classes('w-full flex mt-6').style('height: 75vh'):

		with ui.scroll_area().classes('w-1/5 h-full'):
			# table_of_content = ui.table(columns=columns, rows=[], row_key='index')
			table_of_contents = ui.aggrid({
				'defaultColDef': {'flex': 1},
				'columnDefs': [
					{'headerName': 'Title', 'field': 'title'},
				],
				'rowData': [],
				'rowSelection': 'multiple',
			}).classes('max-h-40')

			for i in range(200):
				table_of_contents.call_api_method('updateRowData', {'add': [
													  {'title': 'test'}
				]})
				#table_of_content.add_rows({'index': index, 'title': title})

		with ui.scroll_area().classes('w-0 pl-14 grow h-full'):
			...

			#for i in range(200):
				#ui.label("hello")

	footer()

