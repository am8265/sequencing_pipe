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

public class extraction_metric : base_cycle_metric {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal extraction_metric(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_metricsPINVOKE.extraction_metric_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(extraction_metric obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~extraction_metric() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_extraction_metric(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public extraction_metric() : this(c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_0(), true) {
  }

  public extraction_metric(extraction_metric_header header) : this(c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_1(extraction_metric_header.getCPtr(header)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public extraction_metric(uint lane, uint tile, uint cycle, ulong date_time, ushort_vector intensities_p90, float_vector focus_scores) : this(c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_2(lane, tile, cycle, date_time, ushort_vector.getCPtr(intensities_p90), float_vector.getCPtr(focus_scores)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  static private global::System.IntPtr SwigConstructextraction_metric(uint lane, uint tile, uint cycle, ulong date_time, ushort[] intensities_p90, float[] focus_scores, uint channel_count) {
    unsafe{ fixed ( ushort* swig_ptrTo_intensities_p90 = intensities_p90 ) {
    unsafe{ fixed ( float* swig_ptrTo_focus_scores = focus_scores ) {
    return c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_3(lane, tile, cycle, date_time, (global::System.IntPtr)swig_ptrTo_intensities_p90, (global::System.IntPtr)swig_ptrTo_focus_scores, channel_count);
    }}
    }}
  }

  public extraction_metric(uint lane, uint tile, uint cycle, ulong date_time, ushort[] intensities_p90, float[] focus_scores, uint channel_count) : this(extraction_metric.SwigConstructextraction_metric(lane, tile, cycle, date_time, intensities_p90, focus_scores, channel_count), true) {
  }

  static private global::System.IntPtr SwigConstructextraction_metric(uint lane, uint tile, uint cycle, ulong date_time, ushort[] intensities_p90, float[] focus_scores) {
    unsafe{ fixed ( ushort* swig_ptrTo_intensities_p90 = intensities_p90 ) {
    unsafe{ fixed ( float* swig_ptrTo_focus_scores = focus_scores ) {
    return c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_4(lane, tile, cycle, date_time, (global::System.IntPtr)swig_ptrTo_intensities_p90, (global::System.IntPtr)swig_ptrTo_focus_scores);
    }}
    }}
  }

  public extraction_metric(uint lane, uint tile, uint cycle, ulong date_time, ushort[] intensities_p90, float[] focus_scores) : this(extraction_metric.SwigConstructextraction_metric(lane, tile, cycle, date_time, intensities_p90, focus_scores), true) {
  }

  public extraction_metric(uint lane, uint tile, uint cycle, csharp_date_time date_time, ushort_vector intensities_p90, float_vector focus_scores) : this(c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_5(lane, tile, cycle, csharp_date_time.getCPtr(date_time), ushort_vector.getCPtr(intensities_p90), float_vector.getCPtr(focus_scores)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  static private global::System.IntPtr SwigConstructextraction_metric(uint lane, uint tile, uint cycle, csharp_date_time date_time, ushort[] intensities_p90, float[] focus_scores, uint channel_count) {
    unsafe{ fixed ( ushort* swig_ptrTo_intensities_p90 = intensities_p90 ) {
    unsafe{ fixed ( float* swig_ptrTo_focus_scores = focus_scores ) {
    return c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_6(lane, tile, cycle, csharp_date_time.getCPtr(date_time), (global::System.IntPtr)swig_ptrTo_intensities_p90, (global::System.IntPtr)swig_ptrTo_focus_scores, channel_count);
    }}
    }}
  }

  public extraction_metric(uint lane, uint tile, uint cycle, csharp_date_time date_time, ushort[] intensities_p90, float[] focus_scores, uint channel_count) : this(extraction_metric.SwigConstructextraction_metric(lane, tile, cycle, date_time, intensities_p90, focus_scores, channel_count), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  static private global::System.IntPtr SwigConstructextraction_metric(uint lane, uint tile, uint cycle, csharp_date_time date_time, ushort[] intensities_p90, float[] focus_scores) {
    unsafe{ fixed ( ushort* swig_ptrTo_intensities_p90 = intensities_p90 ) {
    unsafe{ fixed ( float* swig_ptrTo_focus_scores = focus_scores ) {
    return c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_7(lane, tile, cycle, csharp_date_time.getCPtr(date_time), (global::System.IntPtr)swig_ptrTo_intensities_p90, (global::System.IntPtr)swig_ptrTo_focus_scores);
    }}
    }}
  }

  public extraction_metric(uint lane, uint tile, uint cycle, csharp_date_time date_time, ushort[] intensities_p90, float[] focus_scores) : this(extraction_metric.SwigConstructextraction_metric(lane, tile, cycle, date_time, intensities_p90, focus_scores), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public extraction_metric(uint lane, uint tile, uint cycle, ushort_vector max_intensity_values, float_vector focus_scores) : this(c_csharp_metricsPINVOKE.new_extraction_metric__SWIG_8(lane, tile, cycle, ushort_vector.getCPtr(max_intensity_values), float_vector.getCPtr(focus_scores)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void set(uint lane, uint tile, uint cycle, ulong date_time, ushort_vector max_intensity_values, float_vector focus_scores) {
    c_csharp_metricsPINVOKE.extraction_metric_set__SWIG_0(swigCPtr, lane, tile, cycle, date_time, ushort_vector.getCPtr(max_intensity_values), float_vector.getCPtr(focus_scores));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void set(uint lane, uint tile, uint cycle, ushort_vector max_intensity_values, float_vector focus_scores) {
    c_csharp_metricsPINVOKE.extraction_metric_set__SWIG_1(swigCPtr, lane, tile, cycle, ushort_vector.getCPtr(max_intensity_values), float_vector.getCPtr(focus_scores));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public ulong date_time() {
    ulong ret = c_csharp_metricsPINVOKE.extraction_metric_date_time__SWIG_0(swigCPtr);
    return ret;
    }

  public csharp_date_time date_time_csharp() {
    csharp_date_time ret = new csharp_date_time(c_csharp_metricsPINVOKE.extraction_metric_date_time_csharp(swigCPtr), false);
    return ret;
  }

  public ulong date_time_csharp_raw() {
    ulong ret = c_csharp_metricsPINVOKE.extraction_metric_date_time_csharp_raw(swigCPtr);
    return ret;
    }

  public ushort max_intensity(uint channel) {
    ushort ret = c_csharp_metricsPINVOKE.extraction_metric_max_intensity(swigCPtr, channel);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public float focus_score(uint channel) {
    float ret = c_csharp_metricsPINVOKE.extraction_metric_focus_score(swigCPtr, channel);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public ushort_vector max_intensity_values() {
    ushort_vector ret = new ushort_vector(c_csharp_metricsPINVOKE.extraction_metric_max_intensity_values(swigCPtr), false);
    return ret;
  }

  public float_vector focus_scores() {
    float_vector ret = new float_vector(c_csharp_metricsPINVOKE.extraction_metric_focus_scores(swigCPtr), false);
    return ret;
  }

  public uint channel_count() {
    uint ret = c_csharp_metricsPINVOKE.extraction_metric_channel_count(swigCPtr);
    return ret;
  }

  public void trim(uint channel_count) {
    c_csharp_metricsPINVOKE.extraction_metric_trim(swigCPtr, channel_count);
  }

  public float focusScore(uint channel) {
    float ret = c_csharp_metricsPINVOKE.extraction_metric_focusScore(swigCPtr, channel);
    return ret;
  }

  public ulong dateTime() {
    ulong ret = c_csharp_metricsPINVOKE.extraction_metric_dateTime(swigCPtr);
    return ret;
    }

  public float_vector focusScores() {
    float_vector ret = new float_vector(c_csharp_metricsPINVOKE.extraction_metric_focusScores(swigCPtr), false);
    return ret;
  }

  public void date_time(ulong time) {
    c_csharp_metricsPINVOKE.extraction_metric_date_time__SWIG_1(swigCPtr, time);
  }

  public bool is_any_p90_zero() {
    bool ret = c_csharp_metricsPINVOKE.extraction_metric_is_any_p90_zero(swigCPtr);
    return ret;
  }

  public static string prefix() {
    string ret = c_csharp_metricsPINVOKE.extraction_metric_prefix();
    return ret;
  }

  public static readonly int MAX_CHANNELS = c_csharp_metricsPINVOKE.extraction_metric_MAX_CHANNELS_get();
  public static readonly int TYPE = c_csharp_metricsPINVOKE.extraction_metric_TYPE_get();
  public static readonly int LATEST_VERSION = c_csharp_metricsPINVOKE.extraction_metric_LATEST_VERSION_get();

}

}