FRAMEWORK
---------

Build ``SHARED`` library as Framework Bundle on the OS X and iOS.

If a ``SHARED`` library target has this property set to ``TRUE`` it will be
built as a framework when built on the OS X and iOS.  It will have the
directory structure required for a framework and will be suitable to
be used with the ``-framework`` option

To customize ``Info.plist`` file in the framework, use
:prop_tgt:`MACOSX_FRAMEWORK_INFO_PLIST` target property.

For OS X see also the :prop_tgt:`FRAMEWORK_VERSION` target property.

Example of creation ``dynamicFramework``:

.. code-block:: cmake

  add_library(dynamicFramework SHARED
              dynamicFramework.c
              dynamicFramework.h
  )
  set_target_properties(dynamicFramework PROPERTIES
    FRAMEWORK TRUE
    FRAMEWORK_VERSION C
    MACOSX_FRAMEWORK_IDENTIFIER com.cmake.dynamicFramework
    MACOSX_FRAMEWORK_INFO_PLIST Info.plist
    PUBLIC_HEADER dynamicFramework.h
    XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "iPhone Developer"
  )
