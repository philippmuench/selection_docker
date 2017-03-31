//CHHS treespace.cpp
/* treespace.c
   collection of tree perturbation routines
*/

#include "defines.h"
#include "global.h"
#include "functions.h"

void
ReRootTree (int newroot)
{
  /* reroot tree at newroot.  oldroot forgotten
     The order of branches is not changed.
     Branch lengths, and other parameters for branches are updated.
     Note that node inode needs to be updated if com.oldconP[inode] == 0.
   */
  int oldroot = tree.root, a, b;	/* a->b becomes b->a */

  if (newroot == oldroot)
    {
      return;
    }

  for (b = newroot, a = nodes[b].father; b != oldroot;
       b = a, a = nodes[b].father)
    {
      tree.branches[nodes[b].ibranch][0] = b;
      tree.branches[nodes[b].ibranch][1] = a;
#if (BASEML || CODEML)

      if (a >= com.ns /* && com.method==1 */ )
	{
	  com.oldconP[a] = 0;	/* update the node */
	}

#endif
    }

  tree.root = newroot;
  BranchToNode ();

  for (b = oldroot, a = nodes[b].father; b != newroot;
       b = a, a = nodes[b].father)
    {
      nodes[b].branch = nodes[a].branch;
      nodes[b].label = nodes[a].label;
    }

  nodes[newroot].branch = -1;
  nodes[newroot].label = -1;

#if (CODEML)

  /* omega's are moved in updateconP for NSbranchsites models */
  if (com.model && com.NSsites == 0)
    {
      for (b = oldroot, a = nodes[b].father; b != newroot;
	   b = a, a = nodes[b].father)
	{
	  nodes[b].omega = nodes[a].omega;
	}

      nodes[newroot].omega = -1;
    }

#endif
}

int
NeighborNNI (int i_tree)
{
  /* get the i_tree'th neighboring tree of tree by the nearest neighbor
     interchange (NNI), each tree has 2*(# internal branches) neighbors.
     works with rooted and unrooted binary trees.

     Gives the ip_th neighbor for interior branch ib.
     Involved branches are a..b, a..c, b..d,
     with a..b to be the internal branch.
     swap c with d, with d to be the ip_th son of b
   */
  int i, a, b, c, d, ib = i_tree / 2, ip = i_tree % 2;

  if (tree.nbranch != com.ns * 2 - 2 - (nodes[tree.root].nson == 3))
    {
      error2 ((char *) "err NeighborNNI: multificating tree.");	//CHHS cast to char * to enforce match of datatypes with function call
    }

  /* locate a,b,c,d */
  for (i = 0, a = 0; i < tree.nbranch; i++)
    if (tree.branches[i][1] >= com.ns && a++ == ib)
      {
	break;
      }

  ib = i;
  a = tree.branches[ib][0];
  b = tree.branches[ib][1];
  c = nodes[a].sons[0];

  if (c == b)
    {
      c = nodes[a].sons[1];
    }

  d = nodes[b].sons[ip];

  /* swap nodes c and d */
  tree.branches[nodes[c].ibranch][1] = d;
  tree.branches[nodes[d].ibranch][1] = c;
  BranchToNode ();
  return (0);
}
