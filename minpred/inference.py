"""Small, dependency-flexible TFLite inference wrapper.

Author: Naveen Duhan
"""

from __future__ import annotations

import numpy as np

from minpred.model_store import get_model


def _interpreter_class():
    try:
        from ai_edge_litert.interpreter import Interpreter
        return Interpreter
    except ImportError:
        pass
    try:
        from tflite_runtime.interpreter import Interpreter
        return Interpreter
    except ImportError:
        pass
    try:
        import tensorflow as tf
        return tf.lite.Interpreter
    except ImportError as exc:
        raise RuntimeError(
            "No TFLite runtime is installed. Install MINpred with "
            "`pip install '.[tensorflow]'` or install ai-edge-litert."
        ) from exc


def predict(model_key: str, samples: np.ndarray) -> np.ndarray:
    """Run fixed-batch MINpred TFLite models one sequence at a time."""
    samples = np.asarray(samples)
    if samples.ndim == 2:
        samples = samples[:, np.newaxis, :, np.newaxis]
    if samples.ndim != 4:
        raise ValueError(f"Expected a 2D or 4D feature array; received {samples.shape}")
    if samples.shape[0] == 0:
        return np.empty((0, 0), dtype=np.float32)

    interpreter = _interpreter_class()(model_path=str(get_model(model_key)))
    interpreter.allocate_tensors()
    input_info = interpreter.get_input_details()[0]
    output_info = interpreter.get_output_details()[0]
    expected = tuple(int(value) for value in input_info["shape"][1:])
    if tuple(samples.shape[1:]) != expected:
        raise ValueError(
            f"Model {model_key} expects sample shape {expected}, "
            f"but received {tuple(samples.shape[1:])}."
        )

    outputs = []
    for sample in samples:
        batch = sample[np.newaxis, ...].astype(input_info["dtype"], copy=False)
        interpreter.set_tensor(input_info["index"], batch)
        interpreter.invoke()
        outputs.append(interpreter.get_tensor(output_info["index"])[0].copy())
    return np.asarray(outputs)
