import os
import sys

import click
from flask.cli import with_appcontext


def ptipython_shell(ctx, banner):
    from ptpython.ipython import embed
    from ptpython.repl import run_config
    from ptpython.entry_points.run_ptpython import get_config_and_history_file, create_parser
    a = create_parser().parse_args(args=[])
    config_file, history_file = get_config_and_history_file(a)

    def configure(repl):
        if os.path.exists(config_file):
            run_config(repl, config_file)
    embed(banner1=banner, user_ns=ctx, history_filename=history_file, configure=configure)


def ptpython_shell(ctx, banner):
    from ptpython.repl import embed, run_config
    from ptpython.entry_points.run_ptpython import get_config_and_history_file, create_parser
    a = create_parser().parse_args(args=[])
    config_file, history_file = get_config_and_history_file(a)

    def configure(repl):
        if os.path.exists(config_file):
            run_config(repl, config_file)
    print(banner)
    embed(globals=ctx, history_filename=history_file, configure=configure)


def ipython_shell(ctx, banner):
    from IPython import start_ipython
    from IPython.terminal.ipapp import load_default_config
    config = load_default_config()
    config.TerminalInteractiveShell.banner1 = banner
    start_ipython(argv=[], config=config, user_ns=ctx)


def bpython_shell(ctx, banner):
    from bpython import embed
    embed(banner=banner, locals_=ctx)


def plain_shell(ctx, banner):
    import code
    code.interact(banner=banner, local=ctx)


SHELL_TYPE = ['ptipython', 'ptpython', 'ipython', 'bpython', 'plain']
SHELL_MAP = {
    'ptipython': ptipython_shell,
    'ptpython': ptpython_shell,
    'ipython': ipython_shell,
    'bpython': bpython_shell,
    'plain': plain_shell,
}


@click.command('shell', short_help='Runs a shell in the app context.')
@click.option('--use-shell', type=click.Choice(SHELL_TYPE))
@with_appcontext
def shell_command(use_shell):
    """Runs an interactive Python shell in the context of a given
    Flask application.  The application will populate the default
    namespace of this shell according to it's configuration.

    This is useful for executing small snippets of management code
    without having to manually configuring the application.
    """
    from flask.globals import _app_ctx_stack
    app = _app_ctx_stack.top.app
    banner = 'Python %s on %s\nApp: %s [%s]\nInstance: %s' % (
        sys.version,
        sys.platform,
        app.import_name,
        app.env,
        app.instance_path,
    )
    ctx = {}
    startup = os.environ.get('PYTHONSTARTUP')
    if startup and os.path.isfile(startup):
        with open(startup, 'r') as f:
            eval(compile(f.read(), startup, 'exec'), ctx)

    ctx.update(app.make_shell_context())

    if use_shell:
        try:
            SHELL_MAP.get(use_shell)(ctx, banner)
            return
        except ImportError:
            pass

    for key in SHELL_TYPE:
        try:
            SHELL_MAP.get(key)(ctx, banner)
            return
        except ImportError:
            pass
