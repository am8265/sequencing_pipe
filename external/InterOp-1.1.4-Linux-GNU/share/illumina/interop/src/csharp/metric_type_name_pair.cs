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

public class metric_type_name_pair : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal metric_type_name_pair(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(metric_type_name_pair obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~metric_type_name_pair() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_metric_type_name_pair(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public metric_type_name_pair() : this(c_csharp_metricsPINVOKE.new_metric_type_name_pair__SWIG_0(), true) {
  }

  public metric_type_name_pair(metric_type t, string u) : this(c_csharp_metricsPINVOKE.new_metric_type_name_pair__SWIG_1((int)t, u), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public metric_type_name_pair(metric_type_name_pair p) : this(c_csharp_metricsPINVOKE.new_metric_type_name_pair__SWIG_2(metric_type_name_pair.getCPtr(p)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public metric_type first {
    set {
      c_csharp_metricsPINVOKE.metric_type_name_pair_first_set(swigCPtr, (int)value);
    } 
    get {
      metric_type ret = (metric_type)c_csharp_metricsPINVOKE.metric_type_name_pair_first_get(swigCPtr);
      return ret;
    } 
  }

  public string second {
    set {
      c_csharp_metricsPINVOKE.metric_type_name_pair_second_set(swigCPtr, value);
      if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    } 
    get {
      string ret = c_csharp_metricsPINVOKE.metric_type_name_pair_second_get(swigCPtr);
      if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
      return ret;
    } 
  }

}

}
