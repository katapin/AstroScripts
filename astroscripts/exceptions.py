import os


class TaskError(Exception):
    """Exception to raise in custom scripts."""

    def __init__(self, taskname, custom_message=None, filename=''):
        self.taskname = taskname
        if custom_message:
            self.msg = custom_message
        else:
            self.msg = f"Task '{taskname}' caused an error"
            self.msg += f' during processing file {filename}.' if filename else '.'
        self.filename=filename
        super().__init__(self.msg)


class ExternalTaskError(TaskError):
    """Exception caused due to external program to raise in custom scripts."""

    def __init__(self, taskname: str, *, custom_message: str = None,
                 filename: str | os.PathLike = '', caller: str = ''):
        self.taskname = taskname
        if custom_message:
            self.msg = custom_message
        else:
            all_args=locals()
            extra_info = ["{}='{}'".format(x, all_args[x]) for x in ['filename', 'caller'] if all_args[x] ]
            self.msg = f"Task '{taskname}' finished with error"
            self.msg += ' ({}).'.format(', '.join(extra_info)) if extra_info else '.'
        self.filename = os.fspath(filename)
        self.caller=caller
        super().__init__(taskname, self.msg, self.filename)