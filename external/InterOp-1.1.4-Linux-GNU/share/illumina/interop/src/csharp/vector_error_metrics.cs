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

public class vector_error_metrics : global::System.IDisposable, global::System.Collections.IEnumerable
    , global::System.Collections.Generic.IEnumerable<error_metric>
 {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal vector_error_metrics(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(vector_error_metrics obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~vector_error_metrics() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_vector_error_metrics(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public vector_error_metrics(global::System.Collections.ICollection c) : this() {
    if (c == null)
      throw new global::System.ArgumentNullException("c");
    foreach (error_metric element in c) {
      this.Add(element);
    }
  }

  public bool IsFixedSize {
    get {
      return false;
    }
  }

  public bool IsReadOnly {
    get {
      return false;
    }
  }

  public error_metric this[int index]  {
    get {
      return getitem(index);
    }
    set {
      setitem(index, value);
    }
  }

  public int Capacity {
    get {
      return (int)capacity();
    }
    set {
      if (value < size())
        throw new global::System.ArgumentOutOfRangeException("Capacity");
      reserve((uint)value);
    }
  }

  public int Count {
    get {
      return (int)size();
    }
  }

  public bool IsSynchronized {
    get {
      return false;
    }
  }

  public void CopyTo(error_metric[] array)
  {
    CopyTo(0, array, 0, this.Count);
  }

  public void CopyTo(error_metric[] array, int arrayIndex)
  {
    CopyTo(0, array, arrayIndex, this.Count);
  }

  public void CopyTo(int index, error_metric[] array, int arrayIndex, int count)
  {
    if (array == null)
      throw new global::System.ArgumentNullException("array");
    if (index < 0)
      throw new global::System.ArgumentOutOfRangeException("index", "Value is less than zero");
    if (arrayIndex < 0)
      throw new global::System.ArgumentOutOfRangeException("arrayIndex", "Value is less than zero");
    if (count < 0)
      throw new global::System.ArgumentOutOfRangeException("count", "Value is less than zero");
    if (array.Rank > 1)
      throw new global::System.ArgumentException("Multi dimensional array.", "array");
    if (index+count > this.Count || arrayIndex+count > array.Length)
      throw new global::System.ArgumentException("Number of elements to copy is too large.");
    for (int i=0; i<count; i++)
      array.SetValue(getitemcopy(index+i), arrayIndex+i);
  }

  global::System.Collections.Generic.IEnumerator<error_metric> global::System.Collections.Generic.IEnumerable<error_metric>.GetEnumerator() {
    return new vector_error_metricsEnumerator(this);
  }

  global::System.Collections.IEnumerator global::System.Collections.IEnumerable.GetEnumerator() {
    return new vector_error_metricsEnumerator(this);
  }

  public vector_error_metricsEnumerator GetEnumerator() {
    return new vector_error_metricsEnumerator(this);
  }

  // Type-safe enumerator
  /// Note that the IEnumerator documentation requires an InvalidOperationException to be thrown
  /// whenever the collection is modified. This has been done for changes in the size of the
  /// collection but not when one of the elements of the collection is modified as it is a bit
  /// tricky to detect unmanaged code that modifies the collection under our feet.
  public sealed class vector_error_metricsEnumerator : global::System.Collections.IEnumerator
    , global::System.Collections.Generic.IEnumerator<error_metric>
  {
    private vector_error_metrics collectionRef;
    private int currentIndex;
    private object currentObject;
    private int currentSize;

    public vector_error_metricsEnumerator(vector_error_metrics collection) {
      collectionRef = collection;
      currentIndex = -1;
      currentObject = null;
      currentSize = collectionRef.Count;
    }

    // Type-safe iterator Current
    public error_metric Current {
      get {
        if (currentIndex == -1)
          throw new global::System.InvalidOperationException("Enumeration not started.");
        if (currentIndex > currentSize - 1)
          throw new global::System.InvalidOperationException("Enumeration finished.");
        if (currentObject == null)
          throw new global::System.InvalidOperationException("Collection modified.");
        return (error_metric)currentObject;
      }
    }

    // Type-unsafe IEnumerator.Current
    object global::System.Collections.IEnumerator.Current {
      get {
        return Current;
      }
    }

    public bool MoveNext() {
      int size = collectionRef.Count;
      bool moveOkay = (currentIndex+1 < size) && (size == currentSize);
      if (moveOkay) {
        currentIndex++;
        currentObject = collectionRef[currentIndex];
      } else {
        currentObject = null;
      }
      return moveOkay;
    }

    public void Reset() {
      currentIndex = -1;
      currentObject = null;
      if (collectionRef.Count != currentSize) {
        throw new global::System.InvalidOperationException("Collection modified.");
      }
    }

    public void Dispose() {
        currentIndex = -1;
        currentObject = null;
    }
  }

  public void Clear() {
    c_csharp_metricsPINVOKE.vector_error_metrics_Clear(swigCPtr);
  }

  public void Add(error_metric x) {
    c_csharp_metricsPINVOKE.vector_error_metrics_Add(swigCPtr, error_metric.getCPtr(x));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  private uint size() {
    uint ret = c_csharp_metricsPINVOKE.vector_error_metrics_size(swigCPtr);
    return ret;
  }

  private uint capacity() {
    uint ret = c_csharp_metricsPINVOKE.vector_error_metrics_capacity(swigCPtr);
    return ret;
  }

  private void reserve(uint n) {
    c_csharp_metricsPINVOKE.vector_error_metrics_reserve(swigCPtr, n);
  }

  public vector_error_metrics() : this(c_csharp_metricsPINVOKE.new_vector_error_metrics__SWIG_0(), true) {
  }

  public vector_error_metrics(vector_error_metrics other) : this(c_csharp_metricsPINVOKE.new_vector_error_metrics__SWIG_1(vector_error_metrics.getCPtr(other)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public vector_error_metrics(int capacity) : this(c_csharp_metricsPINVOKE.new_vector_error_metrics__SWIG_2(capacity), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  private error_metric getitemcopy(int index) {
    error_metric ret = new error_metric(c_csharp_metricsPINVOKE.vector_error_metrics_getitemcopy(swigCPtr, index), true);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private error_metric getitem(int index) {
    error_metric ret = new error_metric(c_csharp_metricsPINVOKE.vector_error_metrics_getitem(swigCPtr, index), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void setitem(int index, error_metric val) {
    c_csharp_metricsPINVOKE.vector_error_metrics_setitem(swigCPtr, index, error_metric.getCPtr(val));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void AddRange(vector_error_metrics values) {
    c_csharp_metricsPINVOKE.vector_error_metrics_AddRange(swigCPtr, vector_error_metrics.getCPtr(values));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public vector_error_metrics GetRange(int index, int count) {
    global::System.IntPtr cPtr = c_csharp_metricsPINVOKE.vector_error_metrics_GetRange(swigCPtr, index, count);
    vector_error_metrics ret = (cPtr == global::System.IntPtr.Zero) ? null : new vector_error_metrics(cPtr, true);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Insert(int index, error_metric x) {
    c_csharp_metricsPINVOKE.vector_error_metrics_Insert(swigCPtr, index, error_metric.getCPtr(x));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void InsertRange(int index, vector_error_metrics values) {
    c_csharp_metricsPINVOKE.vector_error_metrics_InsertRange(swigCPtr, index, vector_error_metrics.getCPtr(values));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveAt(int index) {
    c_csharp_metricsPINVOKE.vector_error_metrics_RemoveAt(swigCPtr, index);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveRange(int index, int count) {
    c_csharp_metricsPINVOKE.vector_error_metrics_RemoveRange(swigCPtr, index, count);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static vector_error_metrics Repeat(error_metric value, int count) {
    global::System.IntPtr cPtr = c_csharp_metricsPINVOKE.vector_error_metrics_Repeat(error_metric.getCPtr(value), count);
    vector_error_metrics ret = (cPtr == global::System.IntPtr.Zero) ? null : new vector_error_metrics(cPtr, true);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Reverse() {
    c_csharp_metricsPINVOKE.vector_error_metrics_Reverse__SWIG_0(swigCPtr);
  }

  public void Reverse(int index, int count) {
    c_csharp_metricsPINVOKE.vector_error_metrics_Reverse__SWIG_1(swigCPtr, index, count);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public void SetRange(int index, vector_error_metrics values) {
    c_csharp_metricsPINVOKE.vector_error_metrics_SetRange(swigCPtr, index, vector_error_metrics.getCPtr(values));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

}

}