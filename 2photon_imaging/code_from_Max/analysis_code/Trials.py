import numpy as np


class Trial:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def update(self, **kwargs):
        self.__dict__.update(kwargs)

    def verify(self, **kwargs):
        for key, value in kwargs.items():
            if self.__dict__.get(key, None) != value:
                return False
        return True

    def get(self, key: str):
        if self.__dict__.get(key, None) is not None:
            return self.__dict__.get(key)

    def append(self, func_dict: dict):
        self.__dict__.update({key: func_dict[key](self) for key in func_dict.keys()})


class TrialsCollection:
    def __init__(self, src: list[Trial] = None):
        self.src_data = src if src is not None else []

    def add(self, new_trial: Trial):
        self.src_data.append(new_trial)

    def filter(self, **kwargs):
        result_session = TrialsCollection()
        for trial in self.src_data:
            if trial.verify(**kwargs):
                result_session.add(trial)
        return result_session

    def extract(self, key: str, **kwargs) -> list | np.ndarray:
        result_extraction = []
        for trial in self.src_data:
            if trial.verify(**kwargs):
                result_extraction.append(trial.get(key))
        if len(result_extraction) > 0 and isinstance(result_extraction[0], (np.ndarray, float)):
            result_extraction = np.array(result_extraction)
        return result_extraction

    def extract_super(self, keys: list[str], **kwargs):
        result_dict = {}
        for key in keys:
            result_dict[key] = self.extract(key, **kwargs)
        return result_dict

    def update(self, **kwargs):
        for trial in self.src_data:
            trial.update(**kwargs)

    def __len__(self):
        return len(self.src_data)

    def __getitem__(self, item):
        return self.src_data.__getitem__(item)

    def __add__(self, other):
        if not isinstance(other, TrialsCollection):
            assert other == 0
            return self
        result_session = TrialsCollection(self.src_data + other.src_data)
        return result_session

    def __radd__(self, other):
        if not isinstance(other, TrialsCollection):
            assert other == 0
            return self
        result_session = TrialsCollection(other.src_data + self.src_data)
        return result_session

    def append(self, func_dict: dict):
        for trial in self.src_data:
            trial.append(func_dict)

    def elementwise_update(self, other):
        assert isinstance(other, TrialsCollection)
        assert len(self) == len(other)
        for element_id in range(len(self)):
            self.src_data[element_id].update(**other.src_data[element_id].__dict__)
