//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Summary {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Metrics;
using Illumina.InterOp.Run;
using Illumina.InterOp.RunMetrics;

public class lane_summary_vector : global::System.IDisposable, global::System.Collections.IEnumerable
    , global::System.Collections.Generic.IEnumerable<lane_summary>
 {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal lane_summary_vector(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(lane_summary_vector obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~lane_summary_vector() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_summaryPINVOKE.delete_lane_summary_vector(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public lane_summary_vector(global::System.Collections.ICollection c) : this() {
    if (c == null)
      throw new global::System.ArgumentNullException("c");
    foreach (lane_summary element in c) {
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

  public lane_summary this[int index]  {
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

  public void CopyTo(lane_summary[] array)
  {
    CopyTo(0, array, 0, this.Count);
  }

  public void CopyTo(lane_summary[] array, int arrayIndex)
  {
    CopyTo(0, array, arrayIndex, this.Count);
  }

  public void CopyTo(int index, lane_summary[] array, int arrayIndex, int count)
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

  global::System.Collections.Generic.IEnumerator<lane_summary> global::System.Collections.Generic.IEnumerable<lane_summary>.GetEnumerator() {
    return new lane_summary_vectorEnumerator(this);
  }

  global::System.Collections.IEnumerator global::System.Collections.IEnumerable.GetEnumerator() {
    return new lane_summary_vectorEnumerator(this);
  }

  public lane_summary_vectorEnumerator GetEnumerator() {
    return new lane_summary_vectorEnumerator(this);
  }

  // Type-safe enumerator
  /// Note that the IEnumerator documentation requires an InvalidOperationException to be thrown
  /// whenever the collection is modified. This has been done for changes in the size of the
  /// collection but not when one of the elements of the collection is modified as it is a bit
  /// tricky to detect unmanaged code that modifies the collection under our feet.
  public sealed class lane_summary_vectorEnumerator : global::System.Collections.IEnumerator
    , global::System.Collections.Generic.IEnumerator<lane_summary>
  {
    private lane_summary_vector collectionRef;
    private int currentIndex;
    private object currentObject;
    private int currentSize;

    public lane_summary_vectorEnumerator(lane_summary_vector collection) {
      collectionRef = collection;
      currentIndex = -1;
      currentObject = null;
      currentSize = collectionRef.Count;
    }

    // Type-safe iterator Current
    public lane_summary Current {
      get {
        if (currentIndex == -1)
          throw new global::System.InvalidOperationException("Enumeration not started.");
        if (currentIndex > currentSize - 1)
          throw new global::System.InvalidOperationException("Enumeration finished.");
        if (currentObject == null)
          throw new global::System.InvalidOperationException("Collection modified.");
        return (lane_summary)currentObject;
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
    c_csharp_summaryPINVOKE.lane_summary_vector_Clear(swigCPtr);
  }

  public void Add(lane_summary x) {
    c_csharp_summaryPINVOKE.lane_summary_vector_Add(swigCPtr, lane_summary.getCPtr(x));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  private uint size() {
    uint ret = c_csharp_summaryPINVOKE.lane_summary_vector_size(swigCPtr);
    return ret;
  }

  private uint capacity() {
    uint ret = c_csharp_summaryPINVOKE.lane_summary_vector_capacity(swigCPtr);
    return ret;
  }

  private void reserve(uint n) {
    c_csharp_summaryPINVOKE.lane_summary_vector_reserve(swigCPtr, n);
  }

  public lane_summary_vector() : this(c_csharp_summaryPINVOKE.new_lane_summary_vector__SWIG_0(), true) {
  }

  public lane_summary_vector(lane_summary_vector other) : this(c_csharp_summaryPINVOKE.new_lane_summary_vector__SWIG_1(lane_summary_vector.getCPtr(other)), true) {
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public lane_summary_vector(int capacity) : this(c_csharp_summaryPINVOKE.new_lane_summary_vector__SWIG_2(capacity), true) {
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  private lane_summary getitemcopy(int index) {
    lane_summary ret = new lane_summary(c_csharp_summaryPINVOKE.lane_summary_vector_getitemcopy(swigCPtr, index), true);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private lane_summary getitem(int index) {
    lane_summary ret = new lane_summary(c_csharp_summaryPINVOKE.lane_summary_vector_getitem(swigCPtr, index), false);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void setitem(int index, lane_summary val) {
    c_csharp_summaryPINVOKE.lane_summary_vector_setitem(swigCPtr, index, lane_summary.getCPtr(val));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void AddRange(lane_summary_vector values) {
    c_csharp_summaryPINVOKE.lane_summary_vector_AddRange(swigCPtr, lane_summary_vector.getCPtr(values));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public lane_summary_vector GetRange(int index, int count) {
    global::System.IntPtr cPtr = c_csharp_summaryPINVOKE.lane_summary_vector_GetRange(swigCPtr, index, count);
    lane_summary_vector ret = (cPtr == global::System.IntPtr.Zero) ? null : new lane_summary_vector(cPtr, true);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Insert(int index, lane_summary x) {
    c_csharp_summaryPINVOKE.lane_summary_vector_Insert(swigCPtr, index, lane_summary.getCPtr(x));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void InsertRange(int index, lane_summary_vector values) {
    c_csharp_summaryPINVOKE.lane_summary_vector_InsertRange(swigCPtr, index, lane_summary_vector.getCPtr(values));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveAt(int index) {
    c_csharp_summaryPINVOKE.lane_summary_vector_RemoveAt(swigCPtr, index);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveRange(int index, int count) {
    c_csharp_summaryPINVOKE.lane_summary_vector_RemoveRange(swigCPtr, index, count);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public static lane_summary_vector Repeat(lane_summary value, int count) {
    global::System.IntPtr cPtr = c_csharp_summaryPINVOKE.lane_summary_vector_Repeat(lane_summary.getCPtr(value), count);
    lane_summary_vector ret = (cPtr == global::System.IntPtr.Zero) ? null : new lane_summary_vector(cPtr, true);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Reverse() {
    c_csharp_summaryPINVOKE.lane_summary_vector_Reverse__SWIG_0(swigCPtr);
  }

  public void Reverse(int index, int count) {
    c_csharp_summaryPINVOKE.lane_summary_vector_Reverse__SWIG_1(swigCPtr, index, count);
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

  public void SetRange(int index, lane_summary_vector values) {
    c_csharp_summaryPINVOKE.lane_summary_vector_SetRange(swigCPtr, index, lane_summary_vector.getCPtr(values));
    if (c_csharp_summaryPINVOKE.SWIGPendingException.Pending) throw c_csharp_summaryPINVOKE.SWIGPendingException.Retrieve();
  }

}

}