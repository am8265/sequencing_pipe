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

public class base_read_metric : base_metric {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal base_read_metric(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_metricsPINVOKE.base_read_metric_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(base_read_metric obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~base_read_metric() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_base_read_metric(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public base_read_metric(uint lane, uint tile, uint read) : this(c_csharp_metricsPINVOKE.new_base_read_metric(lane, tile, read), true) {
  }

  public new void set_base(uint lane, uint tile) {
    c_csharp_metricsPINVOKE.base_read_metric_set_base__SWIG_0(swigCPtr, lane, tile);
  }

  public void set_base(uint lane, uint tile, uint read) {
    c_csharp_metricsPINVOKE.base_read_metric_set_base__SWIG_1(swigCPtr, lane, tile, read);
  }

  public uint read() {
    uint ret = c_csharp_metricsPINVOKE.base_read_metric_read(swigCPtr);
    return ret;
  }

  public new ulong id() {
    ulong ret = c_csharp_metricsPINVOKE.base_read_metric_id(swigCPtr);
    return ret;
    }

  public new static ulong create_id(ulong lane, ulong tile, ulong read) {
    ulong ret = c_csharp_metricsPINVOKE.base_read_metric_create_id(lane, tile, read);
    return ret;
    }

  public static ulong read_from_id(ulong id) {
    ulong ret = c_csharp_metricsPINVOKE.base_read_metric_read_from_id(id);
    return ret;
    }

}

}
