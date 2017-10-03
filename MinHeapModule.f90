module MinHeapModule
   implicit none
   public

   ! Constants that require recompilation:
   integer, parameter, private :: dp=kind(0.0d0), sp=kind(0.0)
   !integer, parameter, public :: wp = dp ! Double precision
   integer, parameter, public :: wp = sp ! Single precision

   ! Element of the priority queue.
   type EventType
      real (wp) :: time
      integer :: elementIndex = 0
   end type EventType

   type, public :: MinHeap

      type (EventType), allocatable :: priorityQueue (:)
      integer :: heapSize = 0
      integer, allocatable :: positionInHeap (:)
   
   end type MinHeap

contains

   subroutine initializeHeap(heap,size)
      type (MinHeap), intent(inout) :: heap
      integer, intent(in) :: size

      allocate(heap%priorityQueue(size+1)) ! We use size+1 to avoid overflow when referring to the child of some element
      allocate(heap%positionInHeap(size))
      call resetHeap(heap)

   end subroutine initializeHeap

   subroutine insertInHeap(heap,elementIndex,time) ! AD: Made this routine safe to call even if already in heap
      type (MinHeap), intent(inout) :: heap
      integer, intent(in) :: elementIndex
      real (wp), intent(in) :: time

      integer :: parent, child, tempElementIndex
      real (wp) :: tempTime
      
      if(heap%positionInHeap(elementIndex)/=0) then ! This is already in the heap
         call  updateHeap(heap,elementIndex,time)
         return
      end if

      heap%heapSize = heap%heapSize+1
      heap%priorityQueue(heap%heapSize)%time = time
      heap%priorityQueue(heap%heapSize)%elementIndex = elementIndex
      ! QY: Fortran ignores fractional part for integer division so I didn't distinguish between odd and even cases.   
      parent = heap%heapSize/2
      child = heap%heapSize
      heap%positionInHeap(elementIndex) = heap%heapSize
      SiftUpLoop: do
         if (parent==0) exit SiftUploop
         if (heap%priorityQueue(parent)%time < heap%priorityQueue(child)%time) exit SiftUpLoop
         tempTime = heap%priorityQueue(child)%time 
         tempelementIndex = heap%priorityQueue(child)%elementIndex
         heap%priorityQueue(child)%time = heap%priorityQueue(parent)%time
         heap%priorityQueue(child)%elementIndex = heap%priorityQueue(parent)%elementIndex
         heap%priorityQueue(parent)%time = tempTime
         heap%priorityQueue(parent)%elementIndex = tempelementIndex
         heap%positionInHeap(heap%priorityQueue(child)%elementIndex)=child
         heap%positionInHeap(heap%priorityQueue(parent)%elementIndex)=parent
         child = parent
         parent = parent/2
      end do SiftUpLoop
   end subroutine insertInHeap

   subroutine deleteFromHeap(heap,elementIndex) ! AD: Made this routine safe to call on previously deleted entries
      type (MinHeap), intent(inout) :: heap
      integer, intent(in) :: elementIndex

      integer :: tempIndex, parent, self, last
      real (wp) :: tempTime

      if(heap%positionInHeap(elementIndex)==0) return ! Not in queue!
      if (heap%heapSize==0) then
         write (*,*) elementIndex
         write (*,*) heap%priorityQueue(1)%time,heap%priorityQueue(2)%time,heap%priorityQueue(3)%time
         stop "Deleting from an empty heap"
      end if
      tempIndex = heap%positionInHeap(elementIndex)

      heap%positionInHeap(elementIndex)=0
      if (tempIndex/=heap%heapSize) then
         last = heap%priorityQueue(heap%heapSize)%elementIndex
         heap%priorityQueue(tempIndex)%time = heap%priorityQueue(heap%heapSize)%time
         heap%priorityQueue(tempIndex)%elementIndex = last
         heap%positionInHeap(last)=tempIndex
      end if
      heap%priorityQueue(heap%heapSize)%time = huge(1.0_wp)
      heap%priorityQueue(heap%heapSize)%elementIndex = 0
      heap%heapSize = heap%heapSize - 1

      if ((heap%heapSize>0).and.(tempIndex<=heap%heapSize)) then
         if(heap%priorityQueue(tempIndex)%elementIndex==0) then
            write (*,*) "Heap size is:", heap%heapSize, "cell index is:", elementIndex, "original last one is:", last 
            stop "Deleting heap value of an element that is not in the heap"
         end if
         call modifyHeap(heap,tempIndex)
         if(heap%priorityQueue(tempIndex)%elementIndex==0) stop "Problem arising after modifyHeap"
      end if
   end subroutine deleteFromHeap


   subroutine updateHeap(heap,elementIndex,time)
      type (MinHeap), intent(inout) :: heap
      integer, intent(in) :: elementIndex
      real (wp), intent(in) :: time

      integer :: tempIndex

      tempIndex = heap%positionInHeap(elementIndex)
      if(tempIndex==0) stop "Trying to update heap value that is not in the heap"
      ! QY: we replace the time by a new one. There's no need to replace the elementIndex since it's the same element.
      heap%priorityQueue(tempIndex)%time = time
      call modifyHeap(heap,tempIndex)

   end subroutine updateHeap


   ! QY: new subroutine, used in delete and renew
   subroutine modifyHeap(heap,tempIndex)
      type (MinHeap), intent(inout) :: heap
      integer, intent(in) :: tempIndex

      integer :: parent, leftChild, rightChild, self, minimum
      integer :: tempelementIndex, elementIndex
      real (wp) :: tempTime

      parent = tempIndex/2
      self = tempIndex
      leftChild = tempIndex*2
      tempTime = heap%priorityQueue(tempIndex)%time
      elementIndex = heap%priorityQueue(tempIndex)%elementIndex
      if(elementIndex==0) stop "Trying to modify value that is not actually in the heap"
      
      tempelementIndex = elementIndex
      SiftUpLoop: do
         if (parent==0) exit SiftUpLoop
         if (heap%priorityQueue(parent)%time < tempTime) exit SiftUpLoop
         heap%priorityQueue(self)%time = heap%priorityQueue(parent)%time
         heap%priorityQueue(self)%elementIndex = heap%priorityQueue(parent)%elementIndex
         heap%priorityQueue(parent)%time = tempTime
         heap%priorityQueue(parent)%elementIndex = tempelementIndex
         heap%positionInHeap(heap%priorityQueue(parent)%elementIndex)=parent
         heap%positionInHeap(heap%priorityQueue(self)%elementIndex)=self
         self = parent
         parent = parent/2
      end do SiftUpLoop

      ! In case the node is not sifted up at all, then we may need to sift it down.
      if (self == tempIndex) then
         SiftDownLoop: do
            if (leftChild > heap%heapSize) exit SiftDownLoop
            rightChild=leftChild+1
            minimum = self
            if (heap%priorityQueue(leftChild)%time < tempTime) minimum = leftChild
            if (heap%priorityQueue(rightChild)%time < heap%priorityQueue(minimum)%time) minimum = rightChild
            if (minimum == self) exit SiftDownLoop
            if (minimum == leftChild) then
               tempelementIndex = heap%priorityQueue(self)%elementIndex
               heap%priorityQueue(self)%time = heap%priorityQueue(leftChild)%time
               heap%priorityQueue(self)%elementIndex = heap%priorityQueue(leftChild)%elementIndex
               heap%priorityQueue(leftChild)%time = tempTime
               heap%priorityQueue(leftChild)%elementIndex = tempelementIndex
               heap%positionInHeap(heap%priorityQueue(leftChild)%elementIndex)=leftChild
               heap%positionInHeap(heap%priorityQueue(self)%elementIndex)=self
               self = leftChild
               leftChild = self*2
            else
               tempelementIndex = heap%priorityQueue(self)%elementIndex
               heap%priorityQueue(self)%time = heap%priorityQueue(rightChild)%time
               heap%priorityQueue(self)%elementIndex = heap%priorityQueue(rightChild)%elementIndex
               heap%priorityQueue(rightChild)%time = tempTime
               heap%priorityQueue(rightChild)%elementIndex = tempelementIndex
               heap%positionInHeap(heap%priorityQueue(rightChild)%elementIndex)=rightChild
               heap%positionInHeap(heap%priorityQueue(self)%elementIndex)=self
               self = rightChild
               leftChild = self*2
            end if

         end do SiftDownLoop
      end if
      heap%positionInHeap(elementIndex)=self

   end subroutine modifyHeap


   subroutine resetHeap(heap)
      type (MinHeap), intent(inout) :: heap

      ! In our speciifc case this is the beginning of the time step

      heap%priorityQueue(:)%time = huge(1.0_wp) ! Times of events for each element relative to the last time the queue was reset (time origin)
      heap%priorityQueue(:)%elementIndex = 0
      heap%heapSize = 0
      heap%positionInHeap = 0
      
   end subroutine resetHeap


   subroutine testHeap(heap)
      type (MinHeap), intent(inout) :: heap

      integer :: i, leftChild, rightChild

      do i = 1, heap%heapSize
         leftChild = i*2
         rightChild = i*2+1
         if (leftChild<=heap%heapSize) then
            if (heap%priorityQueue(leftChild)%time < heap%priorityQueue(i)%time) stop 'Wrong heap'
         end if
         if (rightChild<=heap%heapSize) then
            if (heap%priorityQueue(rightChild)%time < heap%priorityQueue(i)%time) stop 'Wrong heap'
         end if
      end do
      do i = 1, heap%heapSize
         if (heap%positionInHeap(heap%priorityQueue(i)%elementIndex)/=i) stop 'Problem in heap'
      end do
   end subroutine testHeap


end module MinHeapModule
