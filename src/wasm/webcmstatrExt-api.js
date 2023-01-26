
const validateSampleSize = (elm) => {
  let isValid = false;
  const iValue = parseInt(elm.value);
  if(iValue != parseFloat(elm.value)) {
    elm.setCustomValidity("Must be a positive integer");
  } else if(iValue < 3) {
    elm.setCustomValidity("Must be at least 3");
  } else {
    elm.setCustomValidity("");
    isValid = true;
  }
  elm.reportValidity();
  return isValid ? iValue : "invalid";
};

const validateSignificance = (elm) => {
  let isValid = false;
  const fValue = parseFloat(elm.value);
  if(fValue < 0.0001 || fValue > 0.5) {
    elm.setCustomValidity("Must be between 0.0001 and 0.5");
  } else {
    elm.setCustomValidity("");
    isValid = true;
  }
  elm.reportValidity();
  return isValid ? fValue : "invalid";
};


const k_equiv_two_sample = (Module, n, m, alpha) => {
  const factor_offset = Module._malloc(8 * 2);
  const factor_array = Module.HEAPF64.subarray(factor_offset/8, factor_offset/8 + 2);
  
  const result = Module.ccall(
      "k_equiv_two_sample", // name of C function
      "number", // return type
      ["number", "number", "number", "number"], // argument types
      [n, m, alpha, factor_offset] // arguments
  );
  
  const k1 = factor_array[0];
  const k2 = factor_array[1];
  
  Module._free(factor_offset);
  
  if(result < 0) {
    throw "Invalid return value when calculating factors.";
  }
  
  return {
    k1: k1,
    k2: k2
  }
};


const power_mean = (Module, n, m, k1, k2) => {
  const steps = 10;
  const mu_1 = 0.;
  const mu_2 = -2.;
  
  let power_x = [];
  let power_y = [];
  
  const mu_delta = (mu_2 - mu_1) / (steps - 1);
  const rejection_rate_offset = Module._malloc(8 * steps);
  const rejection_rate_array = Module.HEAPF64.subarray(rejection_rate_offset / 8,
                                                       rejection_rate_offset / 8 + steps);
  const mu_equiv_offset = Module._malloc(8 * steps);
  const mu_equiv_array = Module.HEAPF64.subarray(mu_equiv_offset / 8,
                                                 mu_equiv_offset / 8 + steps);
  for (let i = 0; i < steps; i++) {
    mu_equiv_array[i] = mu_1 + i * mu_delta;
  }
  
  const result = Module.ccall(
    "power_mean",
    "number",  // return type
    ["number", "number", "number", "number", "number", "number", "number"], // arg types
    [n,        n,        k1,       k2,       steps,    mu_equiv_offset, rejection_rate_offset]
  );
  
  for (let i = 0; i < steps; i++) {
    power_x.push(-mu_equiv_array[i]);
    power_y.push(rejection_rate_array[i]);
  }
  Module._free(mu_equiv_offset);
  Module._free(rejection_rate_offset);
  
  return {
    power_x: power_x,
    power_y: power_y
  };
};

const do_tests = async () => {
  const result = Module.ccall(
    "do_tests",
    "number",  // return type
    [], // arg types
    []
  );
  return result;
};