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

public class base_image_metrics : image_metric_header {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;

  internal base_image_metrics(global::System.IntPtr cPtr, bool cMemoryOwn) : base(c_csharp_metricsPINVOKE.base_image_metrics_SWIGUpcast(cPtr), cMemoryOwn) {
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(base_image_metrics obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~base_image_metrics() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_base_image_metrics(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public base_image_metrics(short version) : this(c_csharp_metricsPINVOKE.new_base_image_metrics__SWIG_0(version), true) {
  }

  public base_image_metrics(image_metric_header header, short version) : this(c_csharp_metricsPINVOKE.new_base_image_metrics__SWIG_1(image_metric_header.getCPtr(header), version), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public base_image_metrics(image_metric_header header) : this(c_csharp_metricsPINVOKE.new_base_image_metrics__SWIG_2(image_metric_header.getCPtr(header)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public base_image_metrics() : this(c_csharp_metricsPINVOKE.new_base_image_metrics__SWIG_3(), true) {
  }

  public base_image_metrics(vector_image_metrics vec, short version, image_metric_header header) : this(c_csharp_metricsPINVOKE.new_base_image_metrics__SWIG_4(vector_image_metrics.getCPtr(vec), version, image_metric_header.getCPtr(header)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public bool data_source_exists() {
    bool ret = c_csharp_metricsPINVOKE.base_image_metrics_data_source_exists__SWIG_0(swigCPtr);
    return ret;
  }

  public void data_source_exists(bool exists) {
    c_csharp_metricsPINVOKE.base_image_metrics_data_source_exists__SWIG_1(swigCPtr, exists);
  }

  public void rebuild_index(bool update_ids) {
    c_csharp_metricsPINVOKE.base_image_metrics_rebuild_index__SWIG_0(swigCPtr, update_ids);
  }

  public void rebuild_index() {
    c_csharp_metricsPINVOKE.base_image_metrics_rebuild_index__SWIG_1(swigCPtr);
  }

  public void resize(uint n) {
    c_csharp_metricsPINVOKE.base_image_metrics_resize(swigCPtr, n);
  }

  public void reserve(uint n) {
    c_csharp_metricsPINVOKE.base_image_metrics_reserve(swigCPtr, n);
  }

  public void trim(uint n) {
    c_csharp_metricsPINVOKE.base_image_metrics_trim(swigCPtr, n);
  }

  public void insert(image_metric metric) {
    c_csharp_metricsPINVOKE.base_image_metrics_insert__SWIG_0(swigCPtr, image_metric.getCPtr(metric));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void insert(ulong id, image_metric metric) {
    c_csharp_metricsPINVOKE.base_image_metrics_insert__SWIG_1(swigCPtr, id, image_metric.getCPtr(metric));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public image_metric at(uint n) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_at(swigCPtr, n), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void set_version(short version) {
    c_csharp_metricsPINVOKE.base_image_metrics_set_version(swigCPtr, version);
  }

  public ulong_vector keys() {
    ulong_vector ret = new ulong_vector(c_csharp_metricsPINVOKE.base_image_metrics_keys(swigCPtr), true);
    return ret;
  }

  public uint_vector lanes() {
    uint_vector ret = new uint_vector(c_csharp_metricsPINVOKE.base_image_metrics_lanes(swigCPtr), true);
    return ret;
  }

  public uint lane_count() {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_lane_count(swigCPtr);
    return ret;
  }

  public uint max_lane() {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_max_lane(swigCPtr);
    return ret;
  }

  public uint_vector tile_numbers_for_lane(uint lane) {
    uint_vector ret = new uint_vector(c_csharp_metricsPINVOKE.base_image_metrics_tile_numbers_for_lane(swigCPtr, lane), true);
    return ret;
  }

  public uint_vector tile_numbers() {
    uint_vector ret = new uint_vector(c_csharp_metricsPINVOKE.base_image_metrics_tile_numbers(swigCPtr), true);
    return ret;
  }

  public vector_image_metrics metrics_for_lane(uint lane) {
    vector_image_metrics ret = new vector_image_metrics(c_csharp_metricsPINVOKE.base_image_metrics_metrics_for_lane__SWIG_0(swigCPtr, lane), true);
    return ret;
  }

  public void metrics_for_lane(vector_image_metrics lane_metrics, uint lane) {
    c_csharp_metricsPINVOKE.base_image_metrics_metrics_for_lane__SWIG_1(swigCPtr, vector_image_metrics.getCPtr(lane_metrics), lane);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public uint_vector cycles() {
    uint_vector ret = new uint_vector(c_csharp_metricsPINVOKE.base_image_metrics_cycles(swigCPtr), true);
    return ret;
  }

  public vector_image_metrics metrics_for_cycle(uint cycle) {
    vector_image_metrics ret = new vector_image_metrics(c_csharp_metricsPINVOKE.base_image_metrics_metrics_for_cycle(swigCPtr, cycle), true);
    return ret;
  }

  public uint size() {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_size(swigCPtr);
    return ret;
  }

  public bool empty() {
    bool ret = c_csharp_metricsPINVOKE.base_image_metrics_empty(swigCPtr);
    return ret;
  }

  public short version() {
    short ret = c_csharp_metricsPINVOKE.base_image_metrics_version(swigCPtr);
    return ret;
  }

  public new void clear() {
    c_csharp_metricsPINVOKE.base_image_metrics_clear(swigCPtr);
  }

  public vector_image_metrics metrics() {
    vector_image_metrics ret = new vector_image_metrics(c_csharp_metricsPINVOKE.base_image_metrics_metrics(swigCPtr), false);
    return ret;
  }

  public static string prefix() {
    string ret = c_csharp_metricsPINVOKE.base_image_metrics_prefix();
    return ret;
  }

  public static string suffix() {
    string ret = c_csharp_metricsPINVOKE.base_image_metrics_suffix();
    return ret;
  }

  public image_metric get_metric_ref(uint lane, uint tile, uint cycle) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric_ref__SWIG_0(swigCPtr, lane, tile, cycle), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public image_metric get_metric_ref(uint lane, uint tile) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric_ref__SWIG_1(swigCPtr, lane, tile), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public image_metric get_metric_ref(ulong key) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric_ref__SWIG_2(swigCPtr, key), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public uint find(uint lane, uint tile, uint cycle) {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_find__SWIG_0(swigCPtr, lane, tile, cycle);
    return ret;
  }

  public uint find(uint lane, uint tile) {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_find__SWIG_1(swigCPtr, lane, tile);
    return ret;
  }

  public uint find(ulong id) {
    uint ret = c_csharp_metricsPINVOKE.base_image_metrics_find__SWIG_2(swigCPtr, id);
    return ret;
  }

  public image_metric get_metric(uint lane, uint tile, uint cycle) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric__SWIG_0(swigCPtr, lane, tile, cycle), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public image_metric get_metric(uint lane, uint tile) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric__SWIG_1(swigCPtr, lane, tile), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public bool has_metric(uint lane, uint tile, uint cycle) {
    bool ret = c_csharp_metricsPINVOKE.base_image_metrics_has_metric__SWIG_0(swigCPtr, lane, tile, cycle);
    return ret;
  }

  public bool has_metric(uint lane, uint tile) {
    bool ret = c_csharp_metricsPINVOKE.base_image_metrics_has_metric__SWIG_1(swigCPtr, lane, tile);
    return ret;
  }

  public image_metric get_metric(ulong key) {
    image_metric ret = new image_metric(c_csharp_metricsPINVOKE.base_image_metrics_get_metric__SWIG_2(swigCPtr, key), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public bool has_metric(ulong id) {
    bool ret = c_csharp_metricsPINVOKE.base_image_metrics_has_metric__SWIG_2(swigCPtr, id);
    return ret;
  }

  public void clear_lookup() {
    c_csharp_metricsPINVOKE.base_image_metrics_clear_lookup(swigCPtr);
  }

  public static readonly int TYPE = c_csharp_metricsPINVOKE.base_image_metrics_TYPE_get();
  public static readonly int LATEST_VERSION = c_csharp_metricsPINVOKE.base_image_metrics_LATEST_VERSION_get();

}

}