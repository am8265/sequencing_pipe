//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Table {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Metrics;
using Illumina.InterOp.Run;
using Illumina.InterOp.RunMetrics;

public class imaging_column : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal imaging_column(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(imaging_column obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~imaging_column() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_tablePINVOKE.delete_imaging_column(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public imaging_column() : this(c_csharp_tablePINVOKE.new_imaging_column__SWIG_0(), true) {
  }

  public imaging_column(uint index, uint offset) : this(c_csharp_tablePINVOKE.new_imaging_column__SWIG_1(index, offset), true) {
  }

  public imaging_column(uint index, uint offset, string_vector sub_columns) : this(c_csharp_tablePINVOKE.new_imaging_column__SWIG_2(index, offset, string_vector.getCPtr(sub_columns)), true) {
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
  }

  public column_id id() {
    column_id ret = (column_id)c_csharp_tablePINVOKE.imaging_column_id__SWIG_0(swigCPtr);
    return ret;
  }

  public string name() {
    string ret = c_csharp_tablePINVOKE.imaging_column_name(swigCPtr);
    return ret;
  }

  public bool has_children() {
    bool ret = c_csharp_tablePINVOKE.imaging_column_has_children(swigCPtr);
    return ret;
  }

  public uint offset() {
    uint ret = c_csharp_tablePINVOKE.imaging_column_offset__SWIG_0(swigCPtr);
    return ret;
  }

  public string_vector subcolumns() {
    string_vector ret = new string_vector(c_csharp_tablePINVOKE.imaging_column_subcolumns(swigCPtr), false);
    return ret;
  }

  public string full_name(uint sub_index) {
    string ret = c_csharp_tablePINVOKE.imaging_column_full_name(swigCPtr, sub_index);
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void offset(uint off) {
    c_csharp_tablePINVOKE.imaging_column_offset__SWIG_1(swigCPtr, off);
  }

  public void id(column_id val) {
    c_csharp_tablePINVOKE.imaging_column_id__SWIG_1(swigCPtr, (int)val);
  }

  public void parse_header_for_id(string header) {
    c_csharp_tablePINVOKE.imaging_column_parse_header_for_id(swigCPtr, header);
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
  }

  public uint size() {
    uint ret = c_csharp_tablePINVOKE.imaging_column_size(swigCPtr);
    return ret;
  }

  public uint column_count() {
    uint ret = c_csharp_tablePINVOKE.imaging_column_column_count(swigCPtr);
    return ret;
  }

  public static string to_header(column_id id) {
    string ret = c_csharp_tablePINVOKE.imaging_column_to_header__SWIG_0((int)id);
    return ret;
  }

  public static string to_header(string name) {
    string ret = c_csharp_tablePINVOKE.imaging_column_to_header__SWIG_1(name);
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string to_name(string header) {
    string ret = c_csharp_tablePINVOKE.imaging_column_to_name__SWIG_0(header);
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string to_name(imaging_column header) {
    string ret = c_csharp_tablePINVOKE.imaging_column_to_name__SWIG_1(imaging_column.getCPtr(header));
    if (c_csharp_tablePINVOKE.SWIGPendingException.Pending) throw c_csharp_tablePINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

}

}
