from django import template
from blog.models import Comments

register = template.Library()


@register.simple_tag()
def comments():
    return Comments.objects.order_by('-id').filter(status="okunmadı")[:5]

