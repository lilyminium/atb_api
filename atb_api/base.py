from pydantic import BaseModel
import numpy as np

from elementable import Elements


class Model(BaseModel):
    class Config:
        validate_all = True
        arbitrary_types_allowed = True
        underscore_attrs_are_private = True
        validate_assignment = True
        allow_population_by_field_name = True
        extra = "forbid"
        json_encoders = {
            np.ndarray: lambda x: x.tolist(),
            tuple: list,
        }
