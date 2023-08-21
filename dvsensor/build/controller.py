from typing import Optional, Callable
from nicegui import background_tasks
import asyncio
import uuid
import utils


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view
		self.current_job_id: Optional[str] = None
		self.current_analysis_job = None

	def init_new_job(self) -> str:
		job_id = str(uuid.uuid4())
		self.current_job_id = job_id
		return job_id

	def clear_job(self):
		self.current_job_id = None
		self.cancel_analysis_job()

	def page_allowed(self, page_kwargs: dict) -> bool:
		if not page_kwargs.get('job_id') or not self.current_job_id or \
				page_kwargs.get('job_id') != self.current_job_id:
			return False
		return True

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		self.model.load_sequence_record(seq_record)
		job_id = self.init_new_job()
		print(f'new job: {job_id}')
		self.view.open_page('metainf', job_id=job_id)

	def init_analysis_job(self, triplet_settings: dict, regions_settings: dict) -> None:
		self.view.open_page('results', job_id=self.current_job_id)

	def start_analysis_job(self, ui_control_functions: dict[str, Callable]) -> dict[str, Callable]:
		print('starting analysis job...')
		self.current_analysis_job = background_tasks.create(self.analysis_job(ui_control_functions))
		job_control_functions: dict[str, Callable] = {
			'cancel_job': self.cancel_analysis_job
		}

		return job_control_functions

	def cancel_analysis_job(self):
		if self.current_analysis_job:
			print('canceling analysis job!')
			self.current_analysis_job.cancel()
			self.current_analysis_job = None

	async def analysis_job(self, ui_control_functions):
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



