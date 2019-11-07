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

public class bar_point : float_point {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal bar_point(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_plotPINVOKE.bar_point_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(bar_point obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~bar_point() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_plotPINVOKE.delete_bar_point(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public bar_point(float x, float height, float width) : this(c_csharp_plotPINVOKE.new_bar_point__SWIG_0(x, height, width), true) {
  }

  public bar_point(float x, float height) : this(c_csharp_plotPINVOKE.new_bar_point__SWIG_1(x, height), true) {
  }

  public bar_point(float x) : this(c_csharp_plotPINVOKE.new_bar_point__SWIG_2(x), true) {
  }

  public bar_point() : this(c_csharp_plotPINVOKE.new_bar_point__SWIG_3(), true) {
  }

  public void set(float x, float height, float width) {
    c_csharp_plotPINVOKE.bar_point_set__SWIG_0(swigCPtr, x, height, width);
  }

  public new void set(float x, float height) {
    c_csharp_plotPINVOKE.bar_point_set__SWIG_1(swigCPtr, x, height);
  }

  public float width() {
    float ret = c_csharp_plotPINVOKE.bar_point_width(swigCPtr);
    return ret;
  }

  public new float min_value() {
    float ret = c_csharp_plotPINVOKE.bar_point_min_value(swigCPtr);
    return ret;
  }

}

}
