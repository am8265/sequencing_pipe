//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Summary {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Metrics;
using Illumina.InterOp.Run;
using Illumina.InterOp.RunMetrics;

public class cycle_state_summary : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal cycle_state_summary(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(cycle_state_summary obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~cycle_state_summary() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_summaryPINVOKE.delete_cycle_state_summary(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public cycle_state_summary() : this(c_csharp_summaryPINVOKE.new_cycle_state_summary(), true) {
  }

  public cycle_range extracted_cycle_range() {
    cycle_range ret = new cycle_range(c_csharp_summaryPINVOKE.cycle_state_summary_extracted_cycle_range__SWIG_0(swigCPtr), false);
    return ret;
  }

  public cycle_range called_cycle_range() {
    cycle_range ret = new cycle_range(c_csharp_summaryPINVOKE.cycle_state_summary_called_cycle_range__SWIG_0(swigCPtr), false);
    return ret;
  }

  public cycle_range qscored_cycle_range() {
    cycle_range ret = new cycle_range(c_csharp_summaryPINVOKE.cycle_state_summary_qscored_cycle_range__SWIG_0(swigCPtr), false);
    return ret;
  }

  public cycle_range error_cycle_range() {
    cycle_range ret = new cycle_range(c_csharp_summaryPINVOKE.cycle_state_summary_error_cycle_range__SWIG_0(swigCPtr), false);
    return ret;
  }

  public bool empty() {
    bool ret = c_csharp_summaryPINVOKE.cycle_state_summary_empty(swigCPtr);
    return ret;
  }

  public void extracted_cycle_range(cycle_range val) {
    c_csharp_summaryPINVOKE.cycle_state_summary_extracted_cycle_range__SWIG_1(swigCPtr, cycle_range.getCPtr(val));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void called_cycle_range(cycle_range val) {
    c_csharp_summaryPINVOKE.cycle_state_summary_called_cycle_range__SWIG_1(swigCPtr, cycle_range.getCPtr(val));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void qscored_cycle_range(cycle_range val) {
    c_csharp_summaryPINVOKE.cycle_state_summary_qscored_cycle_range__SWIG_1(swigCPtr, cycle_range.getCPtr(val));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void error_cycle_range(cycle_range val) {
    c_csharp_summaryPINVOKE.cycle_state_summary_error_cycle_range__SWIG_1(swigCPtr, cycle_range.getCPtr(val));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

}

}