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

public class filter_options : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal filter_options(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(filter_options obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~filter_options() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_plotPINVOKE.delete_filter_options(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read, uint cycle, uint tile_number, uint swath, uint section, uint subsample) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_0((int)naming_method, lane, channel, (int)arg3, surface, read, cycle, tile_number, swath, section, subsample), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read, uint cycle, uint tile_number, uint swath, uint section) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_1((int)naming_method, lane, channel, (int)arg3, surface, read, cycle, tile_number, swath, section), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read, uint cycle, uint tile_number, uint swath) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_2((int)naming_method, lane, channel, (int)arg3, surface, read, cycle, tile_number, swath), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read, uint cycle, uint tile_number) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_3((int)naming_method, lane, channel, (int)arg3, surface, read, cycle, tile_number), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read, uint cycle) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_4((int)naming_method, lane, channel, (int)arg3, surface, read, cycle), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface, uint read) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_5((int)naming_method, lane, channel, (int)arg3, surface, read), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3, uint surface) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_6((int)naming_method, lane, channel, (int)arg3, surface), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel, dna_bases arg3) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_7((int)naming_method, lane, channel, (int)arg3), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane, short channel) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_8((int)naming_method, lane, channel), true) {
  }

  public filter_options(tile_naming_method naming_method, uint lane) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_9((int)naming_method, lane), true) {
  }

  public filter_options(tile_naming_method naming_method) : this(c_csharp_plotPINVOKE.new_filter_options__SWIG_10((int)naming_method), true) {
  }

  public void reset() {
    c_csharp_plotPINVOKE.filter_options_reset(swigCPtr);
  }

  public void validate(metric_type type, info run_info, bool check_ignored) {
    c_csharp_plotPINVOKE.filter_options_validate__SWIG_0(swigCPtr, (int)type, info.getCPtr(run_info), check_ignored);
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
  }

  public void validate(metric_type type, info run_info) {
    c_csharp_plotPINVOKE.filter_options_validate__SWIG_1(swigCPtr, (int)type, info.getCPtr(run_info));
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
  }

  public bool all_channels(metric_type type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_channels__SWIG_0(swigCPtr, (int)type);
    return ret;
  }

  public bool all_channels() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_channels__SWIG_1(swigCPtr);
    return ret;
  }

  public bool all_bases(metric_type type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_bases__SWIG_0(swigCPtr, (int)type);
    return ret;
  }

  public bool all_bases() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_bases__SWIG_1(swigCPtr);
    return ret;
  }

  public bool all_reads() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_reads(swigCPtr);
    return ret;
  }

  public bool all_cycles() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_cycles(swigCPtr);
    return ret;
  }

  public bool all_lanes() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_lanes(swigCPtr);
    return ret;
  }

  public bool all_tile_numbers() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_tile_numbers(swigCPtr);
    return ret;
  }

  public bool all_swaths() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_swaths(swigCPtr);
    return ret;
  }

  public bool all_sections() {
    bool ret = c_csharp_plotPINVOKE.filter_options_all_sections(swigCPtr);
    return ret;
  }

  public bool is_specific_read(metric_type type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_is_specific_read__SWIG_0(swigCPtr, (int)type);
    return ret;
  }

  public bool is_specific_read() {
    bool ret = c_csharp_plotPINVOKE.filter_options_is_specific_read__SWIG_1(swigCPtr);
    return ret;
  }

  public bool is_specific_surface() {
    bool ret = c_csharp_plotPINVOKE.filter_options_is_specific_surface(swigCPtr);
    return ret;
  }

  public bool valid_lane(uint lane) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_lane(swigCPtr, lane);
    return ret;
  }

  public bool valid_surface(uint surface) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_surface(swigCPtr, surface);
    return ret;
  }

  public bool valid_read(uint read) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_read(swigCPtr, read);
    return ret;
  }

  public bool valid_cycle(uint cycle) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_cycle(swigCPtr, cycle);
    return ret;
  }

  public bool valid_tile_number(uint tile_number) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_tile_number(swigCPtr, tile_number);
    return ret;
  }

  public bool valid_swath(uint swath) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_swath(swigCPtr, swath);
    return ret;
  }

  public bool valid_section(uint section) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_section(swigCPtr, section);
    return ret;
  }

  public bool valid_channel(short channel) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_channel(swigCPtr, channel);
    return ret;
  }

  public bool valid_base(dna_bases arg0) {
    bool ret = c_csharp_plotPINVOKE.filter_options_valid_base(swigCPtr, (int)arg0);
    return ret;
  }

  public void subsample(uint count) {
    c_csharp_plotPINVOKE.filter_options_subsample__SWIG_0(swigCPtr, count);
  }

  public void tile_naming_method(tile_naming_method naming_method) {
    c_csharp_plotPINVOKE.filter_options_tile_naming_method(swigCPtr, (int)naming_method);
  }

  public void channel(short channel) {
    c_csharp_plotPINVOKE.filter_options_channel__SWIG_0(swigCPtr, channel);
  }

  public void dna_base(dna_bases arg0) {
    c_csharp_plotPINVOKE.filter_options_dna_base__SWIG_0(swigCPtr, (int)arg0);
  }

  public void read(uint r) {
    c_csharp_plotPINVOKE.filter_options_read__SWIG_0(swigCPtr, r);
  }

  public void cycle(uint c) {
    c_csharp_plotPINVOKE.filter_options_cycle__SWIG_0(swigCPtr, c);
  }

  public void surface(uint s) {
    c_csharp_plotPINVOKE.filter_options_surface__SWIG_0(swigCPtr, s);
  }

  public void swath(uint s) {
    c_csharp_plotPINVOKE.filter_options_swath(swigCPtr, s);
  }

  public void section(uint s) {
    c_csharp_plotPINVOKE.filter_options_section(swigCPtr, s);
  }

  public void tile_number(uint s) {
    c_csharp_plotPINVOKE.filter_options_tile_number(swigCPtr, s);
  }

  public void lane(uint l) {
    c_csharp_plotPINVOKE.filter_options_lane__SWIG_0(swigCPtr, l);
  }

  public uint lane() {
    uint ret = c_csharp_plotPINVOKE.filter_options_lane__SWIG_1(swigCPtr);
    return ret;
  }

  public short channel() {
    short ret = c_csharp_plotPINVOKE.filter_options_channel__SWIG_1(swigCPtr);
    return ret;
  }

  public dna_bases dna_base() {
    dna_bases ret = (dna_bases)c_csharp_plotPINVOKE.filter_options_dna_base__SWIG_1(swigCPtr);
    return ret;
  }

  public uint read() {
    uint ret = c_csharp_plotPINVOKE.filter_options_read__SWIG_1(swigCPtr);
    return ret;
  }

  public uint cycle() {
    uint ret = c_csharp_plotPINVOKE.filter_options_cycle__SWIG_1(swigCPtr);
    return ret;
  }

  public uint surface() {
    uint ret = c_csharp_plotPINVOKE.filter_options_surface__SWIG_1(swigCPtr);
    return ret;
  }

  public string cycle_description() {
    string ret = c_csharp_plotPINVOKE.filter_options_cycle_description(swigCPtr);
    return ret;
  }

  public string lane_description() {
    string ret = c_csharp_plotPINVOKE.filter_options_lane_description(swigCPtr);
    return ret;
  }

  public string channel_description(string_vector channels) {
    string ret = c_csharp_plotPINVOKE.filter_options_channel_description(swigCPtr, string_vector.getCPtr(channels));
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public string base_description() {
    string ret = c_csharp_plotPINVOKE.filter_options_base_description(swigCPtr);
    return ret;
  }

  public string surface_description() {
    string ret = c_csharp_plotPINVOKE.filter_options_surface_description(swigCPtr);
    return ret;
  }

  public string read_description() {
    string ret = c_csharp_plotPINVOKE.filter_options_read_description(swigCPtr);
    return ret;
  }

  public tile_naming_method naming_method() {
    tile_naming_method ret = (tile_naming_method)c_csharp_plotPINVOKE.filter_options_naming_method(swigCPtr);
    return ret;
  }

  public bool supports_section(plot_types arg0, info info) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_section(swigCPtr, (int)arg0, info.getCPtr(info));
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public bool supports_swath(plot_types arg0) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_swath(swigCPtr, (int)arg0);
    return ret;
  }

  public bool supports_tile(plot_types arg0) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_tile(swigCPtr, (int)arg0);
    return ret;
  }

  public bool supports_all_lanes(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_all_lanes(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_lane(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_lane(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_all_bases(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_all_bases(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_base(metric_type metric_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_base(swigCPtr, (int)metric_type);
    return ret;
  }

  public bool supports_all_channels(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_all_channels(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_channel(metric_type metric_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_channel(swigCPtr, (int)metric_type);
    return ret;
  }

  public bool supports_all_cycles(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_all_cycles(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_cycle(metric_type metric_type, plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_cycle(swigCPtr, (int)metric_type, (int)plot_type);
    return ret;
  }

  public bool supports_all_reads(plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_all_reads(swigCPtr, (int)plot_type);
    return ret;
  }

  public bool supports_read(metric_type metric_type, plot_types plot_type) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_read(swigCPtr, (int)metric_type, (int)plot_type);
    return ret;
  }

  public bool supports_surface(metric_type metric_type, info info) {
    bool ret = c_csharp_plotPINVOKE.filter_options_supports_surface(swigCPtr, (int)metric_type, info.getCPtr(info));
    if (c_csharp_plotPINVOKE.SWIGPendingException.Pending) throw c_csharp_plotPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public uint subsample() {
    uint ret = c_csharp_plotPINVOKE.filter_options_subsample__SWIG_1(swigCPtr);
    return ret;
  }

  public enum UseAll {
    ALL_IDS = 0,
    ALL_CHANNELS = -1,
    ALL_BASES = -1
  }

}

}