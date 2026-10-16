#' fastVARs Color Palette
#'
#' Official color palette for the tidyMacro package, designed for VAR model
#' visualizations including impulse response functions, forecast error variance
#' decompositions, and other time series plots.
#'
#' @format A character vector of 14 hex color codes
#'
#' @details
#' The palette consists of 14 distinct colors optimized for clarity and
#' professional presentation:
#' - Crimson (#910048)
#' - Cyan (#00B0B9)
#' - Blue (#407EC9)
#' - Orange (#FF8200)
#' - Gray (#707372)
#' - Navy (#001E60)
#' - Lime (#78BE20)
#' - Purple (#8031A7)
#' - Olive (#658D1B)
#' - Gold (#F2A900)
#' - Dark Orange (#E35205)
#' - Red (#DA291C)
#' - Black (#231F20)
#' - Light Blue (#009CDE)
#'
#' @format A character vector of 13 hex colour codes.
#'
#' @return A character vector of hex colour codes, usable wherever
#'   \pkg{ggplot2} takes colours.
#'
#' @examples
#' # The first three colours are the ones the package's own plots use.
#' tidyMacro_colors[1:3]
#'
#' \donttest{
#' # Use with ggplot2
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
#'   geom_point() +
#'   scale_color_manual(values = tidyMacro_colors) +
#'   fThemeTidyMacro()
#' }
#'
#' @seealso
#' \code{\link{scale_color_tidyMacro}}, \code{\link{scale_fill_tidyMacro}},
#' \code{\link{fThemeTidyMacro}}
#'
#' @export
tidyMacro_colors <- c(
    '#910048',
    '#E5E5E5',
    '#407EC9',
    '#FF8200',
    '#707372',
    '#00B0B9',
    '#001E60',
    '#78BE20',
    '#8031A7',
    '#658D1B',
    '#F2A900',
    '#E35205',
    '#DA291C',
    '#231F20',
    '#009CDE'
)


#' Color Scale Constructor for tidyMacro Palette
#'
#' Internal function to create discrete color or fill scales using the
#' tidyMacro color palette.
#'
#' @param palette Character vector of hex colors. Default is \code{tidyMacro_colors}.
#' @param discrete Logical. If TRUE, creates a discrete scale. Currently only
#'   discrete scales are supported.
#' @param reverse Logical. If TRUE, reverses the order of colors in the palette.
#'   Default is FALSE.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{discrete_scale}}.
#'
#' @return A ggplot2 scale function.
#'
#' @keywords internal
#' @noRd
tidyMacro_pal <- function(
    palette = tidyMacro_colors,
    discrete = TRUE,
    reverse = FALSE,
    ...
) {
    if (reverse) {
        palette <- rev(palette)
    }

    if (discrete) {
        return(function(n) {
            if (n > length(palette)) {
                warning(
                    "Number of groups exceeds available colors. Colors will be recycled.",
                    call. = FALSE
                )
                return(rep(palette, length.out = n))
            }
            palette[1:n]
        })
    }
}


#' tidyMacro Color Scale for ggplot2
#'
#' Apply the tidyMacro color palette to the color aesthetic in ggplot2 plots.
#' This scale is particularly useful for line plots, impulse response functions,
#' and other visualizations where multiple series need distinct colors.
#'
#' @param palette Character vector of hex colors. Default is \code{tidyMacro_colors}.
#' @param discrete Logical. If TRUE (default), creates a discrete scale.
#' @param reverse Logical. If TRUE, reverses the order of colors. Default is FALSE.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{discrete_scale}}.
#'
#' @return A ggplot2 color scale that can be added to a plot.
#'
#' @details
#' This function creates a discrete color scale using the tidyMacro palette.
#' If the number of groups exceeds the 14 available colors, the palette will
#' be recycled with a warning.
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#'
#' # Use the tidyMacro theme for every plot in the session
#' old <- ggplot2::theme_get()
#' fSetTidyMacroTheme()
#' ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) + geom_point()
#'
#' ggplot2::theme_set(old)   # restore
#' }
#'
#' @seealso
#' \code{\link{scale_fill_tidyMacro}}, \code{\link{tidyMacro_colors}},
#' \code{\link{fThemeTidyMacro}}
#'
#' @importFrom ggplot2 discrete_scale
#'
#' @export
scale_color_tidyMacro <- function(
    palette = tidyMacro_colors,
    discrete = TRUE,
    reverse = FALSE,
    ...
) {
    pal <- tidyMacro_pal(
        palette = palette,
        discrete = discrete,
        reverse = reverse
    )

    if (discrete) {
        discrete_scale("colour", "tidyMacro", palette = pal, ...)
    }
}


#' @rdname scale_color_tidyMacro
#' @export
scale_colour_tidyMacro <- scale_color_tidyMacro


#' tidyMacro Fill Scale for ggplot2
#'
#' Apply the tidyMacro color palette to the fill aesthetic in ggplot2 plots.
#' This scale is particularly useful for bar charts, area plots, and confidence
#' bands in impulse response visualizations.
#'
#' @param palette Character vector of hex colors. Default is \code{tidyMacro_colors}.
#' @param discrete Logical. If TRUE (default), creates a discrete scale.
#' @param reverse Logical. If TRUE, reverses the order of colors. Default is FALSE.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{discrete_scale}}.
#'
#' @return A ggplot2 fill scale that can be added to a plot.
#'
#' @details
#' This function creates a discrete fill scale using the tidyMacro palette.
#' If the number of groups exceeds the 14 available colors, the palette will
#' be recycled with a warning.
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(x = factor(cyl), fill = factor(gear))) +
#'   geom_bar(position = "stack") +
#'   labs(x = "Cylinders", y = "Count", fill = "Gears") +
#'   scale_fill_tidyMacro() +
#'   fThemeTidyMacro()
#'
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(gear))) +
#'   geom_col(position = "dodge") +
#'   scale_fill_tidyMacro(reverse = TRUE) +
#'   fThemeTidyMacro()
#' }
#'
#' @seealso
#' \code{\link{scale_color_tidyMacro}}, \code{\link{tidyMacro_colors}},
#' \code{\link{fThemeTidyMacro}}
#'
#' @importFrom ggplot2 discrete_scale
#'
#' @export
scale_fill_tidyMacro <- function(
    palette = tidyMacro_colors,
    discrete = TRUE,
    reverse = FALSE,
    ...
) {
    pal <- tidyMacro_pal(
        palette = palette,
        discrete = discrete,
        reverse = reverse
    )

    if (discrete) {
        discrete_scale("fill", "tidyMacro", palette = pal, ...)
    }
}


#' Update ggplot2 Theme Defaults to tidyMacro Style
#'
#' Sets the tidyMacro theme and color scales as the default for all subsequent
#' ggplot2 plots in the current R session. This is convenient when creating
#' multiple plots with consistent styling.
#'
#' @param base_size Numeric. Base font size for the theme. Default is 14.
#'
#' @return Invisibly returns NULL. Called for side effects.
#'
#' @details
#' This function modifies the default ggplot2 theme and discrete color/fill
#' scales for the current session. All plots created after calling this function
#' will automatically use the tidyMacro styling unless explicitly overridden.
#'
#' To revert to ggplot2 defaults, use \code{ggplot2::theme_set(ggplot2::theme_gray())}
#' and remove the scale defaults with \code{ggplot2::update_geom_defaults()}.
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#' library(tidyMacro)
#'
#' # Set tidyMacro as default theme for the session
#' fSetTidyMacroTheme()
#'
#' # All subsequent plots will use tidyMacro styling
#' ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
#'   geom_point()
#'
#' # Revert to ggplot2 defaults
#' theme_set(theme_gray())
#' }
#'
#' @seealso
#' \code{\link{fThemeTidyMacro}}, \code{\link{scale_color_tidyMacro}},
#' \code{\link[ggplot2]{theme_set}}
#'
#' @importFrom ggplot2 theme_set update_geom_defaults
#'
#' @export
fSetTidyMacroTheme <- function(base_size = 14) {
    # Set the default theme
    theme_set(fThemeTidyMacro(base_size = base_size))

    # Update default aesthetics for common geoms to use tidyMacro colors
    update_geom_defaults("point", list(colour = tidyMacro_colors[1]))
    update_geom_defaults(
        "line",
        list(colour = tidyMacro_colors[1], linewidth = 0.8)
    )
    update_geom_defaults("bar", list(fill = tidyMacro_colors[1]))
    update_geom_defaults("col", list(fill = tidyMacro_colors[1]))

    message("tidyMacro theme and color scales set as default for this session.")
    invisible(NULL)
}
