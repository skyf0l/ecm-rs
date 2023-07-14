use std::mem::MaybeUninit;

/// ECM parameters.
pub struct EcmParams {
    raw: ecm_sys::__ecm_param_struct,
}

impl EcmParams {
    /// Wraps a raw `__ecm_param_struct`.
    ///
    /// # Safety
    ///
    /// The `__ecm_param_struct` must be initialized.
    pub unsafe fn wrap(raw: ecm_sys::__ecm_param_struct) -> Self {
        Self { raw }
    }

    /// Returns a reference to the raw `__ecm_param_struct`.
    pub fn raw(&self) -> &ecm_sys::__ecm_param_struct {
        &self.raw
    }

    /// Returns a new `EcmParams` with default values.
    pub fn new() -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            ecm_sys::ecm_init(raw.as_mut_ptr());
            Self {
                raw: raw.assume_init(),
            }
        }
    }
}

impl Drop for EcmParams {
    fn drop(&mut self) {
        unsafe {
            ecm_sys::ecm_clear(&mut self.raw);
        }
    }
}
