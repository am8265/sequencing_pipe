//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Comm {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Metrics;
using Illumina.InterOp.Run;

public class paths : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal paths(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(paths obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~paths() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_commPINVOKE.delete_paths(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public static string rta_config(string run_directory, int version) {
    string ret = c_csharp_commPINVOKE.paths_rta_config__SWIG_0(run_directory, version);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string rta_config(string run_directory) {
    string ret = c_csharp_commPINVOKE.paths_rta_config__SWIG_1(run_directory);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string run_info(string run_directory) {
    string ret = c_csharp_commPINVOKE.paths_run_info__SWIG_0(run_directory);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string run_parameters(string run_directory, bool alternate) {
    string ret = c_csharp_commPINVOKE.paths_run_parameters__SWIG_0(run_directory, alternate);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string run_parameters(string run_directory) {
    string ret = c_csharp_commPINVOKE.paths_run_parameters__SWIG_1(run_directory);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string run_info() {
    string ret = c_csharp_commPINVOKE.paths_run_info__SWIG_1();
    return ret;
  }

  public static string run_parameters(bool alternate) {
    string ret = c_csharp_commPINVOKE.paths_run_parameters__SWIG_2(alternate);
    return ret;
  }

  public static string run_parameters() {
    string ret = c_csharp_commPINVOKE.paths_run_parameters__SWIG_3();
    return ret;
  }

  public static string rta_config(int version) {
    string ret = c_csharp_commPINVOKE.paths_rta_config__SWIG_2(version);
    return ret;
  }

  public static string rta_config() {
    string ret = c_csharp_commPINVOKE.paths_rta_config__SWIG_3();
    return ret;
  }

  public static string interop_filename(string run_directory, string prefix, string suffix, uint cycle, bool use_out) {
    string ret = c_csharp_commPINVOKE.paths_interop_filename__SWIG_2(run_directory, prefix, suffix, cycle, use_out);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static string interop_filename(string run_directory, string prefix, string suffix, uint cycle) {
    string ret = c_csharp_commPINVOKE.paths_interop_filename__SWIG_3(run_directory, prefix, suffix, cycle);
    if (c_csharp_commPINVOKE.SWIGPendingException.Pending) throw c_csharp_commPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public paths() : this(c_csharp_commPINVOKE.new_paths(), true) {
  }

}

}
