//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Metrics {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Run;

public class q_collapsed_header : q_score_header {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal q_collapsed_header(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_metricsPINVOKE.q_collapsed_header_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(q_collapsed_header obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~q_collapsed_header() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_q_collapsed_header(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public q_collapsed_header() : this(c_csharp_metricsPINVOKE.new_q_collapsed_header__SWIG_0(), true) {
  }

  public q_collapsed_header(q_score_bin_vector bins) : this(c_csharp_metricsPINVOKE.new_q_collapsed_header__SWIG_1(q_score_bin_vector.getCPtr(bins)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public new static q_collapsed_header default_header() {
    q_collapsed_header ret = new q_collapsed_header(c_csharp_metricsPINVOKE.q_collapsed_header_default_header(), true);
    return ret;
  }

  public new void clear() {
    c_csharp_metricsPINVOKE.q_collapsed_header_clear(swigCPtr);
  }

}

}
