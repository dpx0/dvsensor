import asyncio
import logging
import shutil
import nicegui as ng
import threading
import subprocess
from typing import Any, Awaitable
from interface import PageRouter, SignalHandler, JobFn


class AppController:

	def __init__(self, page_router: PageRouter, signal_handler: SignalHandler) -> None:
		self.page_router: PageRouter = page_router
		self.page_router.controller = self

		self.signal_handler = signal_handler
		self._job_thread: Awaitable | None = None
		self._stop_job_thread_event: threading.Event | None = None
		self._job_async_task: asyncio.Task | None = None

	def open_page(self, route: str, source: str | None = None,
				  data: dict[str, Any] | None = None) -> None:
		self.page_router.open(route, source, data)

	def on_page_open(self, data: dict[str, Any] | None) -> None:
		if (not data or not data.get('job_started')) and self.is_job_running():
			logging.debug(f"terminating job")
			self.terminate_job()

	def start_job(self, source: str, job_fn: JobFn, job_data: dict[str, Any]) -> None:
		self.terminate_job()

		self._job_async_task = ng.background_tasks.create(self.execute_job_fn(job_fn, job_data))
		logging.debug(f'job started: {job_fn}({job_data})')

		job_data['job_started'] = True
		self.open_page(route="/results", source=source, data=job_data)

	async def execute_job_fn(self, job_fn: JobFn, job_data: dict[str, Any]) -> None:
		self._stop_job_thread_event = threading.Event()
		self._job_thread = asyncio.to_thread(job_fn,
											 job_data,
											 self._stop_job_thread_event,
											 self.signal_handler)
		await self._job_thread
		self.terminate_job()
		logging.debug(f"job ended: {job_fn}({job_data})")

	def terminate_job(self) -> None:
		actions = [
			(self._stop_job_thread_event, 'set'),
			(self._job_async_task, 'cancel'),
			(self.signal_handler, 'clear_signals')
		]
		for obj, method in actions:
			if hasattr(obj, method):
				getattr(obj, method)()
		self._job_thread = None
		self._stop_job_thread_event = None
		self._job_async_task = None

	def is_job_running(self) -> bool:
		return bool(self._job_thread or (self._stop_job_thread_event and
										 not self._stop_job_thread_event.is_set()))

	@staticmethod
	def detect_blastn_installation() -> str | None:
		return shutil.which('blastn')

	@staticmethod
	def detect_blastdbcheck_installation() -> str | None:
		return shutil.which('blastdbcheck')

	@staticmethod
	def check_blast_db(blastdbcheck_binary: str, db_path: str, db_name: str) -> bool:
		logging.debug(f'running blastdbcheck -db {db_name} in {db_path}')
		try:
			out = subprocess.run([blastdbcheck_binary, '-db', db_name], cwd=db_path,
								 text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		except FileNotFoundError:
			logging.error('blastdbcheck executable not found')
			return False
		except OSError as e:
			logging.error(f'blastdbcheck call failed: {e}')
			return False
		logging.debug(out.stdout)
		if not out.returncode == 0:
			logging.warning(f'blast database check failed')
			return False
		return True
