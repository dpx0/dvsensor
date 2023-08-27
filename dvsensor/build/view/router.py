"""
This code was taken from the nicegui example 'Single Page App' by falkoschindler and
slightly modified (changes indicated by (+) ),
see https://github.com/zauberzeug/nicegui/tree/main/examples/single_page_app
"""

from typing import Awaitable, Callable, Dict, Union
from nicegui import background_tasks, ui


class RouterFrame(ui.element, component='router_frame.js'):
    pass


class Router:

    def __init__(self) -> None:
        self.routes: Dict[str, Callable] = {}
        self.content: ui.element = None

    def add(self, path: str):
        def decorator(func: Callable):
            self.routes[path] = func
            return func
        return decorator

    # + keyword arguments can be passed to builder functions, see below
    def open(self, target: Union[Callable, str], **kwargs) -> None:
        if isinstance(target, str):
            path = target
            # + always redirect to start page if requested page does not exist
            builder = self.routes.get(target, self.routes['/'])
        else:
            path = {v: k for k, v in self.routes.items()}[target]
            builder = target

        async def build() -> None:
            with self.content:
                # + added code to remove all non-nicegui content (e.g. the sensor window on the results page)
                await ui.run_javascript(f'''
                    var toRemove = document.getElementsByClassName("non-ng");
                    Array.from(toRemove).forEach((element, i) => element.parentNode.removeChild(element));
                    if (window.location.pathname !== "{path}") {{
                        history.pushState({{page: "{path}"}}, "", "{path}");
                    }}
                ''', respond=False)
                # + builder functions are called with optional keyword arguments
                result = builder(**kwargs)
                if isinstance(result, Awaitable):
                    await result
        self.content.clear()
        background_tasks.create(build())

    def frame(self) -> ui.element:
        self.content = RouterFrame().on('open', lambda e: self.open(e.args))
        return self.content
