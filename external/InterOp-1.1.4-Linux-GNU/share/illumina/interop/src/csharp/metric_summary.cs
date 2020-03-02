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

public class metric_summary : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal metric_summary(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(metric_summary obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~metric_summary() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_summaryPINVOKE.delete_metric_summary(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public metric_summary(uint arg0) : this(c_csharp_summaryPINVOKE.new_metric_summary(arg0), true) {
  }

  public float error_rate() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_error_rate__SWIG_0(swigCPtr);
    return ret;
  }

  public float percent_aligned() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_percent_aligned__SWIG_0(swigCPtr);
    return ret;
  }

  public float first_cycle_intensity() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_first_cycle_intensity__SWIG_0(swigCPtr);
    return ret;
  }

  public float percent_gt_q30() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_percent_gt_q30__SWIG_0(swigCPtr);
    return ret;
  }

  public float yield_g() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_yield_g__SWIG_0(swigCPtr);
    return ret;
  }

  public float projected_yield_g() {
    float ret = c_csharp_summaryPINVOKE.metric_summary_projected_yield_g__SWIG_0(swigCPtr);
    return ret;
  }

  public void first_cycle_intensity(float val) {
    c_csharp_summaryPINVOKE.metric_summary_first_cycle_intensity__SWIG_1(swigCPtr, val);
  }

  public void error_rate(float val) {
    c_csharp_summaryPINVOKE.metric_summary_error_rate__SWIG_1(swigCPtr, val);
  }

  public void percent_aligned(float val) {
    c_csharp_summaryPINVOKE.metric_summary_percent_aligned__SWIG_1(swigCPtr, val);
  }

  public void percent_gt_q30(float val) {
    c_csharp_summaryPINVOKE.metric_summary_percent_gt_q30__SWIG_1(swigCPtr, val);
  }

  public void yield_g(float val) {
    c_csharp_summaryPINVOKE.metric_summary_yield_g__SWIG_1(swigCPtr, val);
  }

  public void projected_yield_g(float val) {
    c_csharp_summaryPINVOKE.metric_summary_projected_yield_g__SWIG_1(swigCPtr, val);
  }

  public void resize(uint arg0) {
    c_csharp_summaryPINVOKE.metric_summary_resize(swigCPtr, arg0);
  }

}

}