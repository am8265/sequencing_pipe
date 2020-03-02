//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Plot {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.RunMetrics;
using Illumina.InterOp.Metrics;
using Illumina.InterOp.Run;

public class heatmap_data : chart_data {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal heatmap_data(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_plotPINVOKE.heatmap_data_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(heatmap_data obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~heatmap_data() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_plotPINVOKE.delete_heatmap_data(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public heatmap_data() : this(c_csharp_plotPINVOKE.new_heatmap_data(), true) {
  }

  public float at(uint idx) {
    float ret = c_csharp_plotPINVOKE.heatmap_data_at__SWIG_0(swigCPtr, idx);
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public float at(uint row, uint col) {
    float ret = c_csharp_plotPINVOKE.heatmap_data_at__SWIG_1(swigCPtr, row, col);
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public uint row_count() {
    uint ret = c_csharp_plotPINVOKE.heatmap_data_row_count(swigCPtr);
    return ret;
  }

  public uint column_count() {
    uint ret = c_csharp_plotPINVOKE.heatmap_data_column_count(swigCPtr);
    return ret;
  }

  public uint length() {
    uint ret = c_csharp_plotPINVOKE.heatmap_data_length(swigCPtr);
    return ret;
  }

  public bool empty() {
    bool ret = c_csharp_plotPINVOKE.heatmap_data_empty(swigCPtr);
    return ret;
  }

  public void set_buffer(float[] data) {
    unsafe{ fixed ( float* swig_ptrTo_data = data ) {
    {
      c_csharp_plotPINVOKE.heatmap_data_set_buffer__SWIG_0(swigCPtr, (global::System.IntPtr)swig_ptrTo_data);
      if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    }
    }}
  }

  public void set_buffer(float[] data, uint rows, uint cols, float default_val) {
    unsafe{ fixed ( float* swig_ptrTo_data = data ) {
    {
      c_csharp_plotPINVOKE.heatmap_data_set_buffer__SWIG_1(swigCPtr, (global::System.IntPtr)swig_ptrTo_data, rows, cols, default_val);
    }
    }}
  }

  public void set_buffer(float[] data, uint rows, uint cols) {
    unsafe{ fixed ( float* swig_ptrTo_data = data ) {
    {
      c_csharp_plotPINVOKE.heatmap_data_set_buffer__SWIG_2(swigCPtr, (global::System.IntPtr)swig_ptrTo_data, rows, cols);
    }
    }}
  }

  public void resize(uint rows, uint cols, float default_val) {
    c_csharp_plotPINVOKE.heatmap_data_resize__SWIG_0(swigCPtr, rows, cols, default_val);
  }

  public void resize(uint rows, uint cols) {
    c_csharp_plotPINVOKE.heatmap_data_resize__SWIG_1(swigCPtr, rows, cols);
  }

  public new void clear() {
    c_csharp_plotPINVOKE.heatmap_data_clear(swigCPtr);
  }

  public uint index_of(uint row, uint col) {
    uint ret = c_csharp_plotPINVOKE.heatmap_data_index_of(swigCPtr, row, col);
    return ret;
  }

}

}