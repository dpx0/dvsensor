from typing import Optional, Callable
from nicegui import background_tasks
import asyncio
import uuid
import utils


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view
		self.task_controller = TaskController(self)

	def page_allowed(self, page_kwargs: dict) -> bool:
		if not page_kwargs.get('task_id') or not self.task_controller.current_task_id or \
				page_kwargs.get('task_id') != self.task_controller.current_task_id:
			return False
		return True

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		self.model.load_sequence_record(seq_record)
		task_id = self.task_controller.init_new_task()
		print(f'new task: {task_id}')
		self.view.open_page('metainf', task_id=task_id)


class TaskController:

	def __init__(self, controller: Controller) -> None:

		self.controller: Controller = controller
		self.current_task_id: Optional[str] = None
		self.current_analysis_task = None

	def init_new_task(self) -> str:
		task_id = str(uuid.uuid4())
		self.current_task_id = task_id
		return task_id

	def clear_task(self):
		self.current_task_id = None
		self.cancel_analysis_task()

	def init_analysis_task(self, settings) -> None:
		self.controller.view.open_page('results', task_id=self.current_task_id)

	def start_analysis_task(self, ui_control_functions: dict[str, Callable]) -> dict[str, Callable]:
		print('starting analysis task...')
		self.current_analysis_task = background_tasks.create(self.analysis_task(ui_control_functions))
		task_control_functions: dict[str, Callable] = {
			'cancel_task': self.cancel_analysis_task
		}

		return task_control_functions

	def cancel_analysis_task(self) -> None:
		if self.current_analysis_task:
			print('canceling analysis task!')
			self.current_analysis_task.cancel()
			self.current_analysis_task = None

	async def analysis_task(self, ui_control_functions) -> None:
		add_rows = ui_control_functions['add_rows']
		update_progress = ui_control_functions['update_progress']

		data = [{'position': new_id,
			 'triplet': new_id,
			 'region': new_id,
			 'range': new_id,
			 'percent_gc': new_id,
			 'n_stop_edits': new_id,
			 'off_targets': new_id}
				for new_id in range(100)]

		progress_step = 1.0 / len(data)

		while True:
			await asyncio.sleep(2)
			if not data:
				break
			add_rows([data.pop()])

			progress = update_progress(progress_step)
			print(progress)
			if round(progress * 100) >= 100:
				print("FINISHED")
				break



